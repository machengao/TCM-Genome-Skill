#!/usr/bin/env python3
"""Enumerate tRF windows and paired anticodon-loop halves from predicted mature tRNAs."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from contextlib import ExitStack
from html import escape
from pathlib import Path
from tempfile import TemporaryDirectory

from mature_trna import RNA, digest, fasta


def short_fragments(seq):
    n = len(seq)
    for size in range(19, min(25, n) + 1):
        for i in range(n - size + 1):
            end = i + size
            if i == 0 and end == n:
                kind = "full_length"
            elif i == 0:
                kind = "tRF-5"
            elif end == n:
                kind = "tRF-3"
            else:
                kind = "i-tRF"
            yield i + 1, end, kind


def region(pos, loops, n, his):
    if pos > n - 3:
        return "CCA"
    if pos == 1 and his:
        return "G-1"
    for name, (left, right) in loops.items():
        if left <= pos <= right:
            return name + "_loop"
    return "other"


def read_mature(folder):
    meta = json.loads((folder/"summary.json").read_text(encoding="utf-8"))
    if meta.get("sequence_source") != "trnascan_prediction":
        raise ValueError("Rebuild the mature package from the original tRNAscan-SE prediction run")
    sequences = fasta(folder/"mature_trna.fasta")
    with (folder/"mature_trna.tsv").open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    ids = {r["id"] for r in rows}
    reconstructed = {r["id"] for r in rows if r["status"] == "reconstructed"}
    if (len(ids) != len(rows) or set(sequences) != reconstructed or meta["loci"] != len(rows)
            or meta["reconstructed"] != len(reconstructed)):
        raise ValueError("Mature locus table, FASTA and summary do not match")
    for r in rows:
        seq = r["sequence"]
        if r["input_evidence"] != "trnascan_prediction" or r["status"] not in {"reconstructed", "review"}:
            raise ValueError(f"Unexpected sequence provenance or status at {r['id']}")
        if r["status"] == "review" and not seq:
            continue
        if (r["body_source"] != "prediction" or not seq or sequences.get(r["id"]) != seq
                or not set(seq) <= RNA or not seq.endswith("CCA") or len(seq) != int(r["length"])):
            raise ValueError(f"Invalid mature sequence for {r['id']}")
    return meta, rows


def read_structure(r):
    seq = r["sequence"]
    loops = json.loads(r["model_loops"])
    for left, right in loops.values():
        if not 1 <= left <= right <= len(seq):
            raise ValueError(f"Invalid loop bounds for {r['id']}")
    a = int(r["anticodon_start"])
    if not 1 <= a <= len(seq)-2 or len(r["anticodon"]) != 3 or seq[a-1:a+2] != r["anticodon"]:
        raise ValueError(f"Anticodon does not match mature coordinates for {r['id']}")
    cuts = json.loads(r["anticodon_cut_after"])
    if cuts:
        left, right = loops["anticodon"]
        if (len(cuts) != 5 or cuts != sorted(set(cuts)) or cuts[1:4] != [a, a+1, a+2]
                or not left <= cuts[0] or not cuts[-1] < len(seq) or cuts[-1] > right):
            raise ValueError(f"Invalid annotated cut sites for {r['id']}")
    return loops, cuts


def build(folder, out):
    meta, rows = read_mature(folder)
    fields = ["sequence_id", "parent_id", "class", "start", "end", "length", "start_region",
              "end_region", "pair_id", "cut_after", "sequence", "evidence"]
    unique = {"trf": {}, "t_half": {}}
    totals = Counter()
    half_unresolved = []
    details = []
    with ExitStack() as stack:
        writers = {}
        for lib in unique:
            f = stack.enter_context((out/(lib+".tsv")).open("w", encoding="utf-8", newline=""))
            writers[lib] = csv.writer(f, delimiter="\t")
            writers[lib].writerow(fields)
        for r in rows:
            if r["status"] != "reconstructed":
                continue
            seq = r["sequence"]
            n = len(seq)
            parent = r["id"]
            loops, cuts = read_structure(r)
            local = Counter()
            pairs = []

            def emit(lib, start, end, kind, pair="", cut=""):
                frag = seq[start-1:end]
                if frag not in unique[lib]:
                    unique[lib][frag] = "rna_" + hashlib.sha256(frag.encode("ascii")).hexdigest()
                his = r["g_minus_one"] in {"encoded_inferred","added_inferred"}
                writers[lib].writerow([unique[lib][frag], parent, kind, start, end, len(frag),
                                       region(start, loops, n, his), region(end, loops, n, his),
                                       pair, cut, frag, "enumerated_candidate"])
                local[kind] += 1
                totals[kind] += 1

            for start, end, kind in short_fragments(seq):
                emit("trf", start, end, kind)
            for cut in cuts:
                pair = f"{parent}_cut_{cut}"
                emit("t_half", 1, cut, "5-half", pair, cut)
                emit("t_half", cut+1, n, "3-half", pair, cut)
                pairs.append(f"<li>{cut}|{cut+1}　<code>{seq[:cut]} | {seq[cut:]}</code></li>")
            if not cuts:
                half_unresolved.append(parent)
            detail = "<details><summary>切点及两侧序列</summary><ol>"+"".join(pairs)+"</ol></details>" if cuts else "切点定位待复核"
            details.append(f"<tr><td>{escape(parent)}<br>{escape(r['label'])}</td><td>{n}</td>"
                           f"<td>{local['tRF-5']} / {local['i-tRF']} / {local['tRF-3']}</td>"
                           f"<td>{len(cuts)} {detail}</td></tr>")
    for lib, lookup in unique.items():
        with (out/(lib+".fasta")).open("w", encoding="ascii") as f:
            for seq, key in lookup.items():
                f.write(f">{key}\n{seq}\n")
    summary = dict(label=meta.get("label",""), input_loci=meta["loci"], mature_loci=meta["reconstructed"],
                   trf_windows=sum(totals[k] for k in ("tRF-5","i-tRF","tRF-3","full_length")),
                   unique_trf=len(unique["trf"]), half_pairs=totals["5-half"],
                   half_fragments=totals["5-half"]+totals["3-half"], unique_half=len(unique["t_half"]),
                   by_class=dict(totals), unresolved_cut_loci=half_unresolved,
                   t_half_rule="after aligned Sprinzl positions 33,34,35,36,37",
                   parent_sha256={p:digest(folder/p) for p in ("summary.json","mature_trna.tsv","mature_trna.fasta")})
    (out/"summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
    title = escape(summary["label"] or "成熟 tRNA")+" · 片段候选库"
    html = f"""<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title><style>body{{font:16px/1.65 system-ui,sans-serif;color:#203f37;margin:3vw}}th,td{{text-align:left;padding:.5rem;border-bottom:1px solid #c6d5cf}}table{{width:100%;border-collapse:collapse}}code{{overflow-wrap:anywhere}}a{{color:#056752}}</style>
<h1>{title}</h1><p>{summary["input_loci"]} 个预测位点中有 {summary["mature_loci"]} 个成熟体进入建库。所有记录均为计算枚举候选。</p>
<h2>短 tRF</h2><p>逐位点穷举 19–25 nt。共 {summary["trf_windows"]:,} 个窗口、{summary["unique_trf"]:,} 种不同序列。</p>
<p><a href="trf.tsv">分类、端点与来源 TSV</a>　<a href="trf.fasta">去重序列 FASTA</a></p>
<h2>t-half</h2><p>在反密码子环结构位点 33–37 的每个碱基之后切一次，共 {summary["half_pairs"]} 对、{summary["half_fragments"]} 条来源记录、{summary["unique_half"]} 种不同序列。</p>
<p><a href="t_half.tsv">成对片段 TSV</a>　<a href="t_half.fasta">去重序列 FASTA</a>　<a href="summary.json">统计与输入校验值</a></p>
<table><thead><tr><th>母体位点</th><th>长度 nt</th><th>5′ / 内部 / 3′ 短窗口</th><th>t-half 对数</th></tr></thead><tbody>
{"".join(details)}</tbody></table></html>"""
    (out/"report.html").write_text(html,encoding="utf-8")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, type=Path, help="Predicted mature tRNA package")
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    if args.output.exists():
        p.error("Use a new output directory")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=args.output.parent) as tmp:
        out = Path(tmp)/"fragments"
        out.mkdir()
        build(args.input, out)
        out.rename(args.output)


if __name__ == "__main__":
    main()
