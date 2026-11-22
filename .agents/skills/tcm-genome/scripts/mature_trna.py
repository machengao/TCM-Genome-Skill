#!/usr/bin/env python3
"""Reconstruct canonical mature tRNA candidates from an existing trna.py run."""
import argparse
import csv
import hashlib
import json
import re
import subprocess
from collections import Counter
from dataclasses import dataclass
from html import escape
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.parse import quote

RNA = set("ACGU")
MODEL = Path("/usr/share/trnascan-se/models/TRNAinf.cm")


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def fasta(path):
    seqs = {}
    key = None
    with path.open(encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                key = line[1:].split()[0]
                if key in seqs:
                    raise ValueError(f"Duplicate sequence {key} in {path}")
                seqs[key] = []
            elif line:
                if key is None:
                    raise ValueError(f"Missing FASTA header in {path}")
                seqs[key].append(line.upper())
    return {key: "".join(parts) for key, parts in seqs.items()}


@dataclass
class Locus:
    row: dict
    coords: tuple
    seq: str

    @property
    def kind(self):
        m = re.search(r"(?:tRNA\d*-|trn)([A-Z][a-z]{2})", self.row["label"], re.I)
        return m[1].title() if m else ""


def load(folder, genome):
    meta = json.loads((folder / "summary.json").read_text(encoding="utf-8-sig"))
    seqs = fasta(folder / "trna.fasta")
    with (folder / "trna.tsv").open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if (len({r["id"] for r in rows}) != len(rows) or
            set(seqs) != {r["id"] for r in rows if r["status"] == "usable"} or
            meta["loci"] != len(rows) or meta["extracted"] != len(seqs)):
        raise ValueError(f"Inconsistent FASTA, locus table or summary in {folder}")
    loci = []
    for r in rows:
        coords = []
        seq = seqs.get(r["id"], "").replace("T", "U")
        if seq or r["status"] == "usable":
            spans = json.loads(r["spans"])
            neg = r["strand"] == "-1"
            if r["strand"] not in {"1", "-1"} or not spans or not seq:
                raise ValueError(f"Missing sequence or coordinates for {r['id']}")
            # GenBank feature.parts is in transcript order, GFF extraction concatenates
            # ascending exons and reverse-complements the whole product on minus.
            if meta["input_type"] != "genbank" or r["evidence"] == "trnascan_prediction":
                spans = sorted(spans, reverse=neg)
            for a, b in spans:
                if not 1 <= a <= b:
                    raise ValueError(f"Invalid interval in {r['id']}")
                coords.extend(range(b, a-1, -1) if neg else range(a, b+1))
            if len(coords) != len(seq) or len(set(coords)) != len(seq):
                raise ValueError(f"Length or exon overlap mismatch in {r['id']}")
            if genome is not None:
                g = genome[r["seq_id"]]
                bases = "".join(g[c-1] for c in coords)
                if neg:
                    bases = bases.translate(str.maketrans("ACGTRYMKSWBDHVN", "TGCAYRKMSWVHDBN"))
                if bases.replace("T", "U") != seq:
                    raise ValueError(f"Reference sequence mismatch in {r['id']}")
        loci.append(Locus(r, tuple(coords), seq))
    return meta, loci


def reference_matches(primary, references):
    by_base = {}
    for i, loc in enumerate(references):
        group = by_base.setdefault((loc.row["seq_id"], loc.row["strand"]), {})
        for pos in loc.coords:
            group.setdefault(pos, []).append(i)
    for loc in primary:
        group = by_base.get((loc.row["seq_id"], loc.row["strand"]), {})
        hits = {i for pos in loc.coords for i in group.get(pos, ())}
        yield [references[i] for i in sorted(hits)]


def edited(loc, edits):
    if not edits:
        return loc
    seq = list(loc.seq)
    for pos, e in edits.items():
        if pos in loc.coords:
            j = loc.coords.index(pos)
            if seq[j] != e["from"]:
                raise ValueError(f"Edit disagrees with reference base at {loc.row['seq_id']}:{pos}")
            seq[j] = e["to"]
    return Locus(loc.row, loc.coords, "".join(seq))


def hypotheses(loc):
    if not set(loc.seq) <= RNA:
        return
    starts = [0]
    if loc.kind == "His" and loc.seq.startswith("G"):
        starts.append(1)
    for start in starts:
        for k in range(4):
            if k and not loc.seq.endswith("CCA"[:k]):
                continue
            body = loc.seq[start:len(loc.seq)-k]
            if body:
                yield start, k, body


def align(bodies, args, work):
    if not bodies:
        return {}, {}
    with (work / "body.fa").open("w", encoding="ascii") as f:
        for body, key in bodies.items():
            f.write(f">{key}\n{body}\n")
    path = work / "alignment.sto"
    cmd = [*args.runner, args.cmalign, "--cyk", "--nonbanded", "-g", "--outformat", "Pfam",
           "-o", "alignment.sto", args.model.as_posix(), "body.fa"]
    subprocess.run(cmd, cwd=work, check=True, stdout=subprocess.DEVNULL)
    aligned = {}
    tracks = {}
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith("#=GC "):
            _, key, value = line.split()
            tracks[key] = value
        elif line.startswith("s"):
            key, value = line.split()
            aligned[key] = value
    if set(aligned) != set(bodies.values()):
        raise ValueError("Alignment omitted input sequences")
    return aligned, tracks


def structure(body, seq, tracks):
    rf, ss = tracks["RF"], tracks["SS_cons"]
    cols = [i for i, c in enumerate(rf) if c != "."]
    # Verified tRNAscan-SE general model: seven outer pairs, final unpaired 73,
    # and *** in RF for the anticodon. Refuse models with a different layout.
    if (len(cols) < 16 or "".join(ss[i] for i in cols[:7]) != "("*7 or
            "".join(ss[i] for i in cols[-8:-1]) != ")"*7 or
            rf.count("*") != 3 or rf[cols[-1]] != "+"):
        raise ValueError("Unsupported covariance model architecture")
    ix = {}
    for i, char in enumerate(seq.upper()):
        if char in RNA:
            ix[i] = len(ix)
    if "".join(seq[i].upper() for i in ix) != body:
        raise ValueError("Alignment sequence mismatch")
    if [ix.get(i) for i in cols[:7]] != list(range(7)):
        return None
    # Some U/C discriminators are assigned to an adjacent insertion column.
    # Require exactly one residue beyond the intact terminal acceptor pair.
    if [ix.get(i) for i in cols[-8:-1]] != list(range(len(body)-8, len(body)-1)):
        return None
    anti = [ix.get(i) for i, c in enumerate(rf) if c == "*"]
    if None in anti or anti != list(range(anti[0], anti[0]+3)):
        return None
    return anti[0]


def loop_annotation(seq, tracks, offset):
    """Map D, anticodon and T loops and canonical cuts into 1-based mature coordinates."""
    rf, ss = tracks["RF"], tracks["SS_cons"]
    consensus = [i for i, c in enumerate(rf) if c != "."]
    runs = list(re.finditer(r"_+", "".join(ss[i] for i in consensus)))
    stars = [i for i, c in enumerate(rf) if c == "*"]
    middle = next((run for run in runs if consensus[run.start()] <= stars[0]
                   and consensus[run.end()-1] >= stars[-1]), None)
    if middle is None or middle in (runs[0], runs[-1]):
        raise ValueError("Cannot locate anticodon loop in covariance model")
    positions = {}
    pos = offset
    for col, base in enumerate(seq.upper()):
        if base in RNA:
            pos += 1
            positions[col] = pos
    loops = {}
    for label, run in (("D", runs[0]), ("anticodon", middle), ("T", runs[-1])):
        left, right = consensus[run.start()], consensus[run.end()-1]
        present = [positions[i] for i in range(left,right+1) if i in positions]
        if present:
            loops[label] = [present[0],present[-1]]
    center = consensus.index(stars[0])
    cut_cols = consensus[center-1:center+4]
    cuts = []
    if (len(cut_cols) == 5 and cut_cols[1:4] == stars
            and all(ss[i] == "_" and i in positions for i in cut_cols)):
        cuts = [positions[i] for i in cut_cols]
    return loops, cuts


def evaluate(loc, bodies, aligned, tracks, compartment, expected):
    results = []
    for lead, k, body in hypotheses(loc):
        anti = structure(body, aligned[bodies[body]], tracks)
        if anti is None or (expected and body[anti:anti+3] != expected):
            continue
        is_his = loc.kind == "His"
        if not loc.kind and body[anti:anti+3] == "GUG":
            continue
        if is_his and (body[anti:anti+3] != "GUG" or (not lead and compartment != "nuclear")):
            continue
        # Encoded G-1 is upstream of model column 1. For nuclear His only,
        # a missing G-1 can be added as a predicted processing event.
        prepend = is_his and not lead
        seq = ("G" if prepend else "") + loc.seq + "CCA"[k:]
        coords = ([None] if prepend else []) + list(loc.coords) + [None]*(3-k)
        results.append((seq, coords, k, lead, prepend, anti + 1 + int(is_his)))
    return results


def write_table(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def downstream_match(genome, loc, tail):
    if genome is None or not tail:
        return ""
    reference = genome[loc.row["seq_id"]]
    last = loc.coords[-1]
    if loc.row["strand"] == "1":
        flank = reference[last:last+len(tail)]
    else:
        flank = reference[max(0,last-1-len(tail)):last-1][::-1].translate(str.maketrans("ACGT","TGCA"))
    if len(flank) != len(tail) or set(flank)-set("ACGT"):
        return ""
    for i, (a,b) in enumerate(zip(flank,tail)):
        if a != b:
            return i
    return len(tail)


def rebuild(args, work, out):
    genome = fasta(args.genome) if args.genome else None
    meta, loci = load(args.input, genome)
    if meta["route"] != "trnascan_prediction":
        raise ValueError("Use a tRNAscan-SE prediction run as the primary input")
    # Older fallback runs can include unexported reference rows before predictions.
    loci = [loc for loc in loci if loc.row["evidence"] == "trnascan_prediction"]
    if not args.annotation and (args.input / "reference_annotation").is_dir():
        args.annotation = args.input / "reference_annotation"
    refs = []
    if args.annotation:
        rm, refs = load(args.annotation, genome)
        if rm["route"] != "database_annotation" or rm["compartment"] != meta["compartment"]:
            raise ValueError("Annotation must be a database annotation run for the same compartment")
        if loci and refs and not ({x.row["seq_id"] for x in loci} & {x.row["seq_id"] for x in refs}):
            raise ValueError("Input and annotation lack shared reference sequence identifiers")
    edits = json.loads(args.edits.read_text(encoding="utf-8-sig")) if args.edits else []
    by_locus = {}
    for e in edits:
        group = by_locus.setdefault((e["locus_id"], e["seq_id"]), {})
        pos = e["genome_position"]
        if (not isinstance(pos, int) or pos < 1 or pos in group or not e["source"].strip()
                or e["from"] not in RNA or e["to"] not in RNA or e["from"] == e["to"]):
            raise ValueError(f"Invalid documented edit: {e}")
        group[pos] = e
    used = set()
    bodies = {}
    pairs = list(zip(loci, reference_matches(loci, refs)))
    ref_counts = Counter(ref.row["id"] for _, matches in pairs for ref in matches)
    for loc in loci:
        changes = by_locus.get((loc.row["id"], loc.row["seq_id"]), {})
        for _, _, body in hypotheses(edited(loc, changes)):
            if body not in bodies:
                bodies[body] = f"s{len(bodies)}"
    aligned, tracks = align(bodies, args, work)
    table, mapping = [], []
    with (out / "mature_trna.fasta").open("w", encoding="ascii") as fasta_out:
        for loc, matches in pairs:
            r = loc.row
            ref = matches[0] if len(matches) == 1 else None
            findings = []
            if len(matches)>1 or (ref and ref_counts[ref.row["id"]] != 1):
                findings.append("ambiguous_overlap")
            if ref:
                if ref.coords != loc.coords:
                    findings.append("exon_coordinates_differ")
                if loc.kind and ref.kind and loc.kind != ref.kind:
                    findings.append("identity_differs")
            comparison = ",".join(findings) or ("concordant" if ref else "no_overlap" if args.annotation else "not_compared")
            record = dict(reference_comparison=comparison, model_loops="", anticodon_cut_after="",id=r["id"], seq_id=r["seq_id"], strand=r["strand"], label=r["label"],
                          predicted_or_input_id=r["annotation_id"], input_spans=r["spans"],
                          input_evidence=r["evidence"], reference_id=ref.row["annotation_id"] if ref else "",
                          reference_label=ref.row["label"] if ref else "", reference_spans=ref.row["spans"] if ref else "",
                          input_length=len(loc.seq), input_rna=loc.seq, status="review", reason="", body_source="", alignment_id="",
                          anticodon="", anticodon_start="", length="", encoded_cca_bases="", added_tail="", downstream_match_nt="",
                          g_minus_one="", editing="unassessed", sequence="")
            table.append(record)
            if not loc.seq or set(loc.seq)-RNA:
                record["reason"] = "unusable_or_ambiguous_input"
                continue
            expected = r["anticodon"].upper().replace("T", "U")
            if not re.fullmatch("[ACGU]{3}", expected):
                expected = ""
            edit_key = r["id"], r["seq_id"]
            changes = by_locus.get(edit_key, {})
            solutions = {}
            item = edited(loc, changes)
            for solution in evaluate(item, bodies, aligned, tracks, meta["compartment"], "" if changes else expected):
                seq, coords, k, lead, prepend, start = solution
                if expected and changes and seq[start-1:start+2] != expected:
                    before = "".join(changes[p]["from"] if p in changes else b
                                     for p,b in zip(coords[start-1:start+2],seq[start-1:start+2]))
                    if before != expected:
                        continue
                solutions.setdefault((seq, tuple(coords), start, k), solution)
            if len(solutions) != 1:
                record["reason"] = "ambiguous_boundary" if solutions else "structure_or_boundary_unresolved"
                continue
            seq, coords, k, lead, prepend, start = next(iter(solutions.values()))
            key = bodies[item.seq[lead:len(item.seq)-k]]
            loops, cuts = loop_annotation(aligned[key], tracks, lead + int(prepend))
            edits_at = {}
            for genome_pos, e in changes.items():
                pos = coords.index(genome_pos)
                if pos >= len(seq)-3:
                    raise ValueError(f"Invalid edit for {r['id']} at {pos+1}")
                edits_at[pos] = e
                used.add((*edit_key, genome_pos))
            record.update(status="reconstructed", reason="", body_source="prediction", model_loops=json.dumps(loops),
                          anticodon_cut_after=json.dumps(cuts), anticodon=seq[start-1:start+2],
                          alignment_id=key,
                          anticodon_start=start, length=len(seq), encoded_cca_bases=k, added_tail="CCA"[k:],
                          downstream_match_nt=downstream_match(genome, item, "CCA"[k:]),
                          g_minus_one="encoded_inferred" if lead else "added_inferred" if prepend else "not_applicable",
                          editing="documented_edits_applied" if edits_at else "unassessed", sequence=seq)
            fasta_out.write(f">{r['id']} {r['seq_id']} predicted_mature_RNA\n{seq}\n")
            for j,(base, coord) in enumerate(zip(seq, coords)):
                role = "CCA" if j>=len(seq)-3 else "G-1" if j==0 and (lead or prepend) else "body"
                e = edits_at.get(j, {})
                mapping.append(dict(id=r["id"], mature_pos=j+1, base=base, seq_id=r["seq_id"],
                                    strand=r["strand"], genome_pos=coord, role=role,
                                    origin="edited" if e else "genomic" if coord else "inferred_addition",
                                    before_edit=e.get("from",""), edit_source=e.get("source","")))
    if len(used) != len(edits):
        raise ValueError("Some supplied edits did not apply to a reconstructed locus")
    fields = "id seq_id strand label predicted_or_input_id input_spans input_evidence reference_id reference_label reference_spans reference_comparison input_length input_rna status reason body_source alignment_id anticodon anticodon_start model_loops anticodon_cut_after length encoded_cca_bases added_tail downstream_match_nt g_minus_one editing sequence".split()
    write_table(out/"mature_trna.tsv", fields, table)
    write_table(out/"base_map.tsv", "id mature_pos base seq_id strand genome_pos role origin before_edit edit_source".split(), mapping)
    inputs = [d/name for d in [args.input, args.annotation] if d for name in ("trna.tsv","trna.fasta","summary.json")]
    inputs += [p for p in (args.genome,args.edits) if p]
    if args.wsl_distro:
        model_hash = subprocess.check_output([*args.runner,"sha256sum",args.model.as_posix()]).split()[0].decode("ascii")
    else:
        model_hash = digest(args.model)
    summary = dict(sequence_source="trnascan_prediction", compartment=meta["compartment"], loci=len(table), comparison_loci=len(refs), reconstructed=sum(bool(r["sequence"]) for r in table),
                   unique_sequences=len({r["sequence"] for r in table if r["sequence"]}),
                   review=dict(Counter(r["reason"] for r in table if r["reason"])),
                   input_sha256={str(p.resolve()):digest(p) for p in inputs},
                   alignment_tool=next((line for line in subprocess.check_output([*args.runner,args.cmalign,"-h"],encoding="utf-8").splitlines() if "INFERNAL" in line),""),
                   model=args.model.as_posix(), model_sha256=model_hash, editing_evidence=str(args.edits or ""), label=args.label)
    (out/"summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
    if bodies:
        (work/"alignment.sto").replace(out/"alignment.sto")
        (work/"body.fa").replace(out/"alignment_inputs.fasta")
    report(out, table, summary)


def report(out, rows, s):
    def h(x):
        return escape(str(x))
    reasons = dict(identity_conflict="预测与注释的氨基酸类型冲突", ambiguous_boundary="多个成熟边界均符合结构",
                   structure_or_boundary_unresolved="受体茎、反密码子或末端未能唯一定位",
                   reference_boundaries_disagree="重叠注释的两端均与输入不一致",
                   ambiguous_reference="输入与注释缺乏一对一对应", unusable_or_ambiguous_input="输入不可用或含模糊碱基")
    names = {"prediction": "扫描预测外显子"}
    title = (s["label"]+" · " if s["label"] else "")+"成熟 tRNA 重建"
    html = [f'<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{h(title)}</title>',
            '<style>body{font:16px/1.6 system-ui,sans-serif;color:#173b33;margin:3vw}h1{font-size:2rem}table{width:100%;border-collapse:collapse}td,th{text-align:left;vertical-align:top;padding:.55rem;border-bottom:1px solid #c6d5cf}th{background:#e8f2eb}code{font:14px/1.5 monospace;overflow-wrap:anywhere}a{color:#056752}small{color:#4c625a}section{overflow-x:auto}</style>',
            f'<h1>{h(title)}</h1><p>{s["loci"]} 个输入位点 · {s["reconstructed"]} 条重建记录 · {s["unique_sequences"]} 种不同序列 · {s["loci"]-s["reconstructed"]} 个待复核</p>',
            '<p>按转录方向展示 5′ → 3′ 的计算重建序列。表中“补”指相对所选位点边界补齐的 CCA。边界外相邻的参考碱基有时也与 CCA 相容，实际末端加工需要转录本证据。RNA 编辑仅在提供位点证据时写入。化学修饰未展开。</p>',
            '<p><a href="mature_trna.fasta">RNA FASTA</a>　<a href="mature_trna.tsv">完整位点表</a>　<a href="base_map.tsv">逐碱基坐标</a>　<a href="summary.json">运行记录与校验值</a></p>']
    totals = Counter(r["body_source"] for r in rows if r["sequence"])
    lengths = [r["length"] for r in rows if r["sequence"]]
    if lengths:
        html.append(f'<p>长度 {min(lengths)}–{max(lengths)} nt。'+
                    '，'.join(f'{names[key]} {num} 条' for key,num in totals.items())+'。</p>')
    if s["comparison_loci"]:
        html.append(f'<p>本批次以 {s["loci"]} 个输入位点为范围，另用同一参考的 {s["comparison_loci"]} 个发布注释位点核对。</p>')
    if not s["editing_evidence"]:
        html.append('<p>本次运行以参考碱基重建，未提供位点级 RNA 编辑证据。</p>')
    html.append('<section><table><thead><tr><th>输入位点 / 链</th><th>输入名称 / 基因组坐标</th><th>加工依据</th><th>nt</th><th>成熟序列或复核原因</th></tr></thead><tbody>')
    for r in rows:
        seq = r["sequence"]
        accession = r["seq_id"]
        src = h(names.get(r["body_source"],"待复核"))
        if seq:
            src += f'<br><small>补 {h(r["added_tail"] or "0 nt")}'
            if r["downstream_match_nt"]:
                src += f' · 下游相容 {h(r["downstream_match_nt"])} nt'
            if r["g_minus_one"] != "not_applicable":
                src += " · G−1 "+("来自参考" if r["g_minus_one"]=="encoded_inferred" else "推断添加")
            src += '</small>'
            if r["input_spans"] != r["reference_spans"] and r["reference_id"]:
                src += f'<br><small>对照 {h(r["reference_id"])} {h(r["reference_spans"])}</small>'
        else:
            src += '<br><small>'+h(r["reference_label"])+'</small>'
        if r["reference_comparison"] not in {"concordant", "not_compared"}:
            src += f'<br><small>注释对照 {h(r["reference_comparison"])}</small>'
        value = '<code>'+h(seq)+'</code>' if seq else h(reasons[r["reason"]])
        ref_link = h(accession)
        if re.fullmatch(r"(?:NC_|[A-Z]{1,4})\d+\.\d+",accession):
            ref_link = f'<a href="https://www.ncbi.nlm.nih.gov/nuccore/{quote(accession)}">{h(accession)}</a>'
        html.append(f'<tr><td>{h(r["id"])}<br>{h(r["strand"])}</td><td>{h(r["label"])}<br><small>{ref_link} {h(r["input_spans"])}</small></td><td>{src}</td><td>{h(r["length"])}</td><td>{value}</td></tr>')
    html.append('</tbody></table></section></html>')
    (out/"report.html").write_text("\n".join(html),encoding="utf-8")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--annotation", type=Path, help="Matching trna.py database_annotation run")
    p.add_argument("--genome", type=Path, help="Reference FASTA for per-base genomic verification")
    p.add_argument("--edits", type=Path, help="JSON list of source-supported substitutions")
    p.add_argument("--model", type=Path, default=MODEL)
    p.add_argument("--cmalign", default="cmalign")
    p.add_argument("--wsl-distro", help="Use cmalign and the model from this Windows WSL distribution")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--label", default="")
    args = p.parse_args()
    args.runner = ["wsl.exe","-d",args.wsl_distro,"--exec"] if args.wsl_distro else []
    if not args.wsl_distro:
        args.model = args.model.resolve()
    if args.output.exists():
        p.error("Use a new output directory")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=args.output.parent) as tmp:
        work = Path(tmp)
        out = work/"results"
        out.mkdir()
        rebuild(args, work, out)
        out.rename(args.output)


if __name__ == "__main__":
    main()
