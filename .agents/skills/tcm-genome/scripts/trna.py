#!/usr/bin/env python3
"""Predict tRNA loci with tRNAscan-SE and archive published features for comparison."""

import argparse
from collections import defaultdict
import csv
import gzip
from hashlib import file_digest
import json
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from urllib.parse import unquote
from zipfile import BadZipFile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import ExactPosition
from Bio.SeqRecord import SeqRecord

from genome_download import acquire, cached_assets
from genome_sources import QueryError


def open_text(path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig")
    return path.open(encoding="utf-8-sig")


def locus(seq_id, location, strand, spans, label, anticodon, status, sequence="", annotation_id=""):
    return dict(seq_id=seq_id, location=location, strand=strand, spans=spans, label=label,
                anticodon=anticodon, status=status, sequence=sequence, annotation_id=annotation_id)


def scan_genbank(path, expected, fasta, compartment):
    rows = []
    count = 0
    with open_text(path) as source, fasta.open("w", encoding="ascii") as output:
        for item in SeqIO.parse(source, "genbank"):
            count += 1
            if expected and item.id != expected:
                raise ValueError(f"GenBank returned {item.id}, expected {expected}")
            source_feature = next((feature for feature in item.features if feature.type == "source"), None)
            organelle = " ".join(source_feature.qualifiers.get("organelle", [])) if source_feature else ""
            if organelle and compartment not in organelle.lower():
                raise ValueError(f"GenBank source indicates {organelle}, expected {compartment}")
            SeqIO.write(item, output, "fasta")
            for feature in item.features:
                if feature.type != "tRNA":
                    continue
                q = feature.qualifiers
                location = feature.location
                parts = location.parts if location is not None else []
                spans = [(int(part.start) + 1, int(part.end)) for part in parts]
                precise = parts and all(
                    type(part.start) is ExactPosition and type(part.end) is ExactPosition and not part.ref
                    and 0 <= int(part.start) < int(part.end) <= len(item)
                    for part in parts
                )
                pseudo = "pseudo" in q or "pseudogene" in q
                usable = precise and location.strand in (1, -1) and not pseudo
                sequence = str(feature.extract(item.seq)).upper() if usable else ""
                coordinate_text = ",".join(f"{start}..{end}" for start, end in spans)
                rows.append(locus(
                    item.id, coordinate_text, location.strand if location else None, spans,
                    next(iter(q.get("product") or q.get("gene") or ["tRNA"])),
                    next(iter(q.get("anticodon", [""]))),
                    "usable" if usable else "pseudogene" if pseudo else "review", sequence,
                    next(iter(q.get("locus_tag") or q.get("gene") or [""])),
                ))
    if not count or (expected and count != 1):
        raise ValueError("Invalid or empty GenBank download")
    return rows


def gff_loci(path):
    rows = {}
    children = defaultdict(list)
    header = False
    sequence_ids = set()
    with open_text(path) as file:
        for number, raw in enumerate(file, 1):
            line = raw.rstrip("\r\n")
            if line == "##FASTA":
                break
            if line.startswith("##gff-version 3"):
                header = True
            if not line or line.startswith("#"):
                continue
            seq_id, _, kind, start, end, _, strand, _, raw_attributes = line.split("\t")
            attrs = {key: unquote(value) for key, separator, value in
                     (field.partition("=") for field in raw_attributes.split(";")) if separator}
            interval = (int(start), int(end))
            if not 1 <= interval[0] <= interval[1]:
                raise ValueError(f"Invalid GFF3 interval on line {number}")
            seq_id = unquote(seq_id)
            sequence_ids.add(seq_id)
            if kind.lower() == "exon":
                for parent in attrs.get("Parent", "").split(","):
                    children[parent].append((seq_id, strand, interval))
            if kind.lower() not in {"trna", "pseudogenic_trna"}:
                continue
            identifier = attrs.get("ID") or f"line-{number}"
            if identifier in rows:
                raise ValueError(f"Repeated GFF3 tRNA ID {identifier}; review multipart representation")
            flags = " ".join(attrs.get(k, "") for k in ("pseudo", "pseudogene", "product", "Note")).lower()
            pseudo = (kind.lower() == "pseudogenic_trna" or "pseudogene" in attrs
                      or "pseudogene" in flags or attrs.get("pseudo") == "true")
            rows[identifier] = locus(
                seq_id, f"{start}..{end}", {"+": 1, "-": -1}.get(strand), [interval],
                attrs.get("product") or attrs.get("Name") or "tRNA", attrs.get("anticodon", ""),
                "pseudogene" if pseudo else "usable" if strand in ("+", "-") else "review",
                annotation_id=attrs.get("ID", ""),
            )
    if not header and not sequence_ids:
        raise ValueError("Annotation is not a recognizable GFF3 file")
    for identifier, row in rows.items():
        if children[identifier]:
            strand = "+" if row["strand"] == 1 else "-"
            intervals = children[identifier]
            lo, hi = row["spans"][0]
            spans = sorted({span for _, _, span in intervals})
            if (any(seq_id != row["seq_id"] or sign != strand for seq_id, sign, _ in intervals)
                    or any(not lo <= start <= end <= hi for start, end in spans)
                    or any(left[1] >= right[0] for left, right in zip(spans, spans[1:]))):
                raise ValueError(f"Inconsistent exon coordinates for {identifier}")
            row["spans"] = spans
    return list(rows.values()), sequence_ids


def fasta_extract(fasta, annotation, excluded):
    rows, annotated_ids = annotation
    by_id = defaultdict(list)
    for row in rows:
        if row["seq_id"] not in excluded:
            by_id[row["seq_id"]].append(row)
    seen = set()
    with open_text(fasta) as source:
        for item in SeqIO.parse(source, "fasta"):
            if item.id in seen:
                raise ValueError(f"Duplicate FASTA ID: {item.id}")
            seen.add(item.id)
            for row in by_id.get(item.id, []):
                if not all(1 <= start <= end <= len(item) for start, end in row["spans"]):
                    raise ValueError(f"GFF3 coordinates out of range for {item.id}")
                if row["status"] == "usable":
                    raw = Seq("".join(str(item.seq[lo - 1:hi]) for lo, hi in row["spans"]))
                    row["sequence"] = str(raw.reverse_complement() if row["strand"] == -1 else raw).upper()
    if not seen or annotated_ids - seen:
        raise ValueError(f"FASTA lacks annotation sequence IDs: {sorted(annotated_ids - seen)[:5]}")
    return [row for row in rows if row["seq_id"] not in excluded]


def prediction_input(fasta, excluded, output):
    with open_text(fasta) as source, output.open("w", encoding="ascii") as target:
        count = SeqIO.write((row for row in SeqIO.parse(source, "fasta") if row.id not in excluded), target, "fasta")
    if not count:
        raise ValueError("No sequences remained for tRNAscan-SE")


def predict(fasta, excluded, compartment, scanner, threads, work, wsl_distro):
    source = work / "prediction.fna"
    gff = work / "prediction.gff"
    prediction_input(fasta, excluded, source)
    mode = "-E" if compartment == "nuclear" else "-O"
    prefix = [scanner, mode, "--thread", str(threads)]
    table = work / "prediction.txt"
    if wsl_distro:
        command = ["wsl.exe", "--distribution", wsl_distro, "--cd", str(work.resolve()),
                   "--exec", *prefix, "-j", gff.name, "-o", table.name, source.name]
    else:
        command = [*prefix, "-j", str(gff), "-o", str(table), str(source)]
    try:
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except FileNotFoundError as exc:
        raise ValueError(
            f"The required prediction executable is unavailable: {command[0]}"
        ) from exc
    sys.stderr.buffer.write(result.stdout)
    result.check_returncode()
    if (b"No tRNAs found." in result.stdout.splitlines()
            and gff.stat().st_size == 0 and table.stat().st_size == 0):
        return [], gff
    return fasta_extract(source, gff_loci(gff), set()), gff


def write_loci(deliverable, rows):
    with (deliverable / "trna.tsv").open("w", newline="", encoding="utf-8") as handle:
        columns = ["id", "seq_id", "annotation_id", "location", "strand", "spans", "label",
                   "anticodon", "status", "evidence"]
        writer = csv.DictWriter(handle, columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for i, row in enumerate(rows, 1):
            row["id"] = f"trna_{i:05d}"
            writer.writerow({**row, "spans": json.dumps(row["spans"])})
    records = [SeqRecord(Seq(row["sequence"]), id=row["id"], description=row["seq_id"])
               for row in rows if row["sequence"]]
    with (deliverable / "trna.fasta").open("w", encoding="ascii") as handle:
        SeqIO.write(records, handle, "fasta")
    return len(records)


def run(args, work):
    gb, fasta, annotation, report = args.genbank, args.fasta, args.gff, None
    input_paths = [path for path in (gb, fasta, annotation) if path]
    if args.accession:
        if args.download_dir:
            gb, fasta, annotation, report = cached_assets(args.accession, args.download_dir)
            input_paths = [path for path in (gb, fasta, annotation, report) if path]
        else:
            gb, fasta, annotation, report = acquire(args.accession, work)
    excluded = set()
    selected = set(args.sequence_id)
    if report:
        labeled = set()
        wanted = set()
        for line in report.read_text(encoding="utf-8").splitlines():
            info = json.loads(line)
            kind = (info.get("assignedMoleculeLocationType") or "").lower()
            if kind in {"mitochondrion", "chloroplast", "plastid"}:
                ids = {info[key] for key in ("genbankAccession", "refseqAccession") if info.get(key)}
                labeled.update(ids)
                if kind == args.compartment:
                    wanted.update(ids)
        if args.compartment == "nuclear":
            excluded = labeled
        else:
            selected = selected & wanted if selected else wanted
            if not selected:
                raise ValueError("NCBI sequence report has no matching component for the requested compartment")
    if gb:
        if selected:
            raise ValueError("Use a separate GenBank record when selecting an individual sequence")
        fasta = work / "genbank.fna"
        if args.accession and args.compartment == "nuclear":
            raise ValueError("Use an assembly accession for a remote nuclear genome")
        rows = scan_genbank(gb, args.accession, fasta, args.compartment)
    else:
        if selected:
            with open_text(fasta) as source:
                ids = {row.id for row in SeqIO.parse(source, "fasta")}
            if not ids & selected or ids & selected & excluded:
                raise ValueError("Requested sequence IDs are absent or outside the requested compartment")
            excluded |= ids - selected
        elif args.accession and args.accession.startswith(("GCA_", "GCF_")) and not report and args.compartment != "nuclear":
            raise ValueError("Missing NCBI sequence report; supply explicit --sequence-id values")
        rows = fasta_extract(fasta, gff_loci(annotation), excluded) if annotation else []
    for row in rows:
        row["evidence"] = "database_annotation"
    predicted, predicted_gff = predict(
        fasta, excluded, args.compartment, args.scanner, args.threads, work, args.wsl_distro
    )
    for row in predicted:
        row["evidence"] = "trnascan_prediction"
    deliverable = work / "results"
    deliverable.mkdir()
    extracted = write_loci(deliverable, predicted)
    fingerprints = {}
    for role, path in (("genbank", gb), ("fasta", None if gb else fasta), ("gff", annotation), ("sequence_report", report)):
        if path:
            with path.open("rb") as handle:
                fingerprints[role] = file_digest(handle, "sha256").hexdigest()
    summary = dict(accession=args.accession, compartment=args.compartment, route="trnascan_prediction", loci=len(predicted), reference_loci=len(rows),
                   extracted=extracted, input_type="genbank" if gb else "fasta+gff" if annotation else "fasta",
                   excluded_ids=sorted(excluded), input_sha256=fingerprints,
                   input_paths=[str(path) for path in input_paths])
    if gb or annotation:
        ref_dir = deliverable / "reference_annotation"
        ref_dir.mkdir()
        ref_meta = dict(summary, route="database_annotation", loci=len(rows), extracted=write_loci(ref_dir, rows))
        (ref_dir / "summary.json").write_text(json.dumps(ref_meta, indent=2), encoding="utf-8")
    summary["scanner"] = args.scanner
    summary["prediction_loci"] = len(predicted)
    if args.wsl_distro:
        summary["wsl_distro"] = args.wsl_distro
    summary["model"] = "-E" if args.compartment == "nuclear" else "-O"
    predicted_gff.replace(deliverable / "prediction.gff")
    (work / "prediction.txt").replace(deliverable / "prediction.txt")
    (deliverable / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    deliverable.rename(args.output)
    print(json.dumps(summary))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--accession")
    inputs.add_argument("--genbank", type=Path)
    inputs.add_argument("--fasta", type=Path)
    parser.add_argument("--gff", type=Path)
    parser.add_argument("--download-dir", type=Path, help="Reuse an accession bundle downloaded by genome_download.py")
    parser.add_argument("--sequence-id", action="append", default=[], help="Restrict a FASTA to these sequence IDs")
    parser.add_argument("--compartment", choices=["nuclear", "chloroplast", "mitochondrion"], required=True)
    parser.add_argument("--scanner", help="Explicit local tRNAscan-SE executable")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True, help="New directory for extracted sequences and evidence")
    args = parser.parse_args()
    args.wsl_distro = os.getenv("TCM_TRNASCAN_WSL_DISTRO") if os.name == "nt" and args.scanner is None else None
    args.scanner = args.scanner or "tRNAscan-SE"
    if args.gff and not args.fasta:
        parser.error("--gff requires --fasta")
    if args.download_dir and not args.accession:
        parser.error("--download-dir requires --accession")
    if args.threads < 1 or args.output.exists():
        parser.error("Use a positive thread count and a new output directory")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    try:
        with TemporaryDirectory(dir=args.output.parent, prefix=".trna-") as temp:
            run(args, Path(temp))
    except (OSError, ValueError, KeyError, TypeError, BadZipFile, QueryError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f"tRNA extraction stopped: {exc}\n")


if __name__ == "__main__":
    main()
