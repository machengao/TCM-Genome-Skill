#!/usr/bin/env python3
"""Download selected reference accessions and tune concurrency from transfer throughput."""

import argparse
from collections import deque
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from datetime import datetime, timezone
import gzip
from hashlib import file_digest
import json
from pathlib import Path
import re
import shutil
import subprocess
from tempfile import TemporaryDirectory
from threading import Lock
import time
from urllib.parse import urlencode, urlsplit
from zipfile import BadZipFile, ZipFile

from genome_sources import GWH, QueryError, get_json


ASSEMBLY = re.compile(r"GC[AF]_\d+\.\d+")
GWH_ASSEMBLY = re.compile(r"GWH[A-Z]+\d+(?:\.\d+)?")
NUCLEOTIDE = re.compile(r"(?!GC[AF]_)[A-Z][A-Z_]*\d+\.\d+")
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ASSET_NAMES = {
    "genbank": {"sequence.gb"},
    "fasta": {"genome.fna", "genome.fna.gz"},
    "gff": {"annotation.gff", "annotation.gff.gz"},
    "sequence_report": {"sequence_report.jsonl"},
}


def validate_accession(accession):
    if not any(pattern.fullmatch(accession) for pattern in (ASSEMBLY, GWH_ASSEMBLY, NUCLEOTIDE)):
        raise ValueError(f"Provide a versioned NCBI nucleotide or assembly accession, or a GWH accession: {accession}")
    return accession


def run_transfer(command, target):
    subprocess.run(command, check=True)


def download(url, path, transfer):
    transfer(
        ["curl", "--fail", "--silent", "--show-error", "--location", "--proto", "=https",
         "--proto-redir", "=https", "--retry", "2", "--retry-all-errors", "--connect-timeout", "15",
         "--speed-limit", "1024", "--speed-time", "120", "--output", str(path), url],
        path,
    )
    return path


def gwh_asset(url, path, transfer):
    if not url:
        return None
    uri = urlsplit(url)
    if uri.scheme != "https" or uri.hostname not in {"download.cncb.ac.cn", "ngdc.cncb.ac.cn"}:
        raise ValueError("Unrecognized GWH download host")
    if uri.path.endswith(".gz"):
        path = path.with_suffix(path.suffix + ".gz")
    return download(url, path, transfer)


def acquire(accession, work, transfer=run_transfer):
    """Return GenBank or paired FASTA and GFF with optional molecule assignments."""
    validate_accession(accession)
    if ASSEMBLY.fullmatch(accession):
        if not shutil.which("datasets"):
            raise ValueError("NCBI Datasets CLI is required for a versioned NCBI assembly accession")
        archive = work / "datasets.zip"
        transfer(
            ["datasets", "download", "genome", "accession", accession,
             "--include", "genome,gff3,seq-report", "--no-progressbar", "--filename", str(archive)],
            archive,
        )
        base = f"ncbi_dataset/data/{accession}/"
        with ZipFile(archive) as package:
            files = {name for name in package.namelist() if name.startswith(base)}
            genomes = [name for name in files if name.endswith("_genomic.fna")]
            if len(genomes) != 1:
                raise ValueError("NCBI package does not contain one FASTA for the exact assembly version")
            assets = [(genomes[0], work / "genome.fna"), (base + "genomic.gff", work / "annotation.gff"),
                      (base + "sequence_report.jsonl", work / "sequence_report.jsonl")]
            for name, target in assets:
                if name in files:
                    with package.open(name) as source, target.open("wb") as output:
                        shutil.copyfileobj(source, output)
        archive.unlink()
        return None, *(target if target.exists() else None for _, target in assets)
    if GWH_ASSEMBLY.fullmatch(accession):
        detail = get_json(f"{GWH}/api/public/assembly/{accession}")
        if detail.get("assemblyAccession") != accession:
            raise ValueError("GWH returned a different assembly accession")
        fasta = gwh_asset(detail.get("ftpPathDna"), work / "genome.fna", transfer)
        if not fasta:
            raise ValueError("GWH did not provide a FASTA URL")
        return None, fasta, gwh_asset(detail.get("ftpPathGff"), work / "annotation.gff", transfer), None
    query = urlencode(dict(db="nuccore", id=accession, rettype="gbwithparts", retmode="text"))
    path = download(f"{EFETCH}?{query}", work / "sequence.gb", transfer)
    return path, None, None, None


def cached_assets(accession, directory):
    """Open a previously completed bundle by its exact accession."""
    folder = directory / validate_accession(accession)
    manifest = json.loads((folder / "manifest.json").read_text(encoding="utf-8"))
    files = manifest.get("files") if isinstance(manifest, dict) else None
    if (not isinstance(files, dict) or manifest.get("accession") != accession
            or not files.keys() <= ASSET_NAMES.keys()
            or ("genbank" in files) == ("fasta" in files)):
        raise ValueError(f"Invalid download manifest for {accession}")
    result = []
    for role, names in ASSET_NAMES.items():
        if role not in files:
            result.append(None)
            continue
        asset = files[role]
        if not isinstance(asset, dict) or asset.get("name") not in names:
            raise ValueError(f"Invalid cached {role} name for {accession}")
        path = folder / asset["name"]
        if path.is_symlink() or not path.is_file() or path.stat().st_size != asset.get("bytes"):
            raise ValueError(f"Cached {role} file is missing or changed for {accession}: {path}")
        result.append(path)
    return tuple(result)


def store_accession(accession, directory, transfer):
    final = directory / accession
    if final.exists():
        cached_assets(accession, directory)
        return "reused"
    with TemporaryDirectory(dir=directory, prefix=f".{accession}-") as scratch:
        staged = Path(scratch) / "reference"
        staged.mkdir()
        assets = acquire(accession, staged, transfer)
        gb, fasta, _, _ = assets
        reference = gb or fasta
        if not reference:
            raise ValueError(f"Missing reference sequence for {accession}")
        expected = b"LOCUS " if gb else b">"
        opener = gzip.open if reference.suffix == ".gz" else open
        with opener(reference, "rb") as source:
            if source.read(len(expected)) != expected:
                raise ValueError(f"Unexpected reference sequence format for {accession}")
        files = {}
        for role, path in zip(ASSET_NAMES, assets):
            if path:
                with path.open("rb") as source:
                    files[role] = {"name": path.name, "bytes": path.stat().st_size,
                                   "sha256": file_digest(source, "sha256").hexdigest()}
        (staged / "manifest.json").write_text(
            json.dumps({"accession": accession, "files": files}, indent=2), encoding="utf-8"
        )
        try:
            staged.rename(final)
        except OSError:
            if not final.is_dir():
                raise
            cached_assets(accession, directory)
            return "reused"
    return "downloaded"


def file_size(path):
    try:
        return path.stat().st_size
    except FileNotFoundError:
        return 0


class TransferMeter:
    def __init__(self):
        self.targets = set()
        self.finished_bytes = 0
        self.lock = Lock()
        self.efetch_lock = Lock()
        self.last_efetch = 0.0

    def run(self, command, target):
        if command[-1].startswith(EFETCH + "?"):
            with self.efetch_lock:
                time.sleep(max(0, 0.35 - (time.monotonic() - self.last_efetch)))
                self.last_efetch = time.monotonic()
        with self.lock:
            self.targets.add(target)
        try:
            result = subprocess.run(command, capture_output=True)
            if result.returncode:
                message = (result.stderr or result.stdout).decode("utf-8", errors="replace").strip()
                raise ValueError(f"{command[0]} exited with code {result.returncode}: {message[-700:]}")
        finally:
            with self.lock:
                self.finished_bytes += file_size(target)
                self.targets.remove(target)

    def sample(self):
        with self.lock:
            return self.finished_bytes + sum(map(file_size, self.targets)), len(self.targets)


class AdaptiveConcurrency:
    def __init__(self, limit):
        self.limit = limit
        self.workers = 1
        self.rate = None
        self.probing = False
        self.cooldown = 0

    def back_off(self):
        self.workers = max(1, self.workers - 1)
        self.probing = False
        self.rate = None
        self.cooldown = 2

    def observe(self, rate, queued):
        if self.probing:
            if rate >= self.rate * 1.15:
                self.rate = rate
                self.probing = False
            else:
                self.back_off()
            return
        if self.workers > 1 and self.rate is not None and rate < self.rate * 0.7:
            self.back_off()
            return
        if rate <= 0:
            return
        self.rate = rate if self.rate is None else (self.rate + rate) / 2
        if self.cooldown:
            self.cooldown -= 1
        elif queued and self.workers < self.limit:
            self.workers += 1
            self.probing = True


def download_batch(accessions, directory, report_path, limit, interval=5.0):
    meter = TransferMeter()
    concurrency = AdaptiveConcurrency(limit)
    pending = deque(accessions)
    running = {}
    report = {"directory": str(directory), "max_downloads": limit, "concurrency": 1,
              "transfer_samples": [], "accessions": {item: {"status": "queued"} for item in accessions}}

    def checkpoint():
        # Replace the report only once the full snapshot has been written.
        with TemporaryDirectory(dir=report_path.parent, prefix=".download-report-") as scratch:
            path = Path(scratch) / "report.json"
            path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
            path.replace(report_path)

    checkpoint()
    before, active = meter.sample()
    sampled = time.monotonic()
    with ThreadPoolExecutor(max_workers=limit) as pool:
        while pending or running:
            launched = False
            while pending and len(running) < concurrency.workers:
                accession = pending.popleft()
                running[pool.submit(store_accession, accession, directory, meter.run)] = accession
                report["accessions"][accession] = {"status": "running"}
                launched = True
            if launched:
                checkpoint()
            done, _ = wait(running, timeout=max(0, sampled + interval - time.monotonic()),
                           return_when=FIRST_COMPLETED)
            changed = bool(done)
            cap = concurrency.workers
            for future in done:
                accession = running.pop(future)
                try:
                    status = {"status": future.result()}
                except (OSError, ValueError, KeyError, TypeError, EOFError, BadZipFile, QueryError) as exc:
                    status = {"status": "error", "error": str(exc)}
                    concurrency.back_off()
                report["accessions"][accession] = status
                print(f"{accession}: {status['status']}", flush=True)
            now = time.monotonic()
            if now - sampled >= interval:
                total, current_active = meter.sample()
                # Require all slots to have active transfers at both measurement points.
                if active == current_active == cap == concurrency.workers:
                    rate = max(0, total - before) / (now - sampled)
                    concurrency.observe(rate, bool(pending))
                    report["transfer_samples"].append(
                        {"downloads": cap, "bytes_per_second": round(rate), "next_limit": concurrency.workers}
                    )
                    del report["transfer_samples"][:-20]
                    changed = True
                before, active, sampled = total, current_active, now
            report["concurrency"] = concurrency.workers
            if changed:
                checkpoint()
    return int(any(row["status"] == "error" for row in report["accessions"].values()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("accessions", nargs="*", help="Selected, versioned reference accessions")
    parser.add_argument("--accession-file", type=Path, help="Text file with one selected accession per line")
    parser.add_argument("--output", type=Path, required=True, help="Directory for reusable reference bundles")
    parser.add_argument("--report", type=Path, help="New JSON progress report path")
    parser.add_argument("--max-downloads", type=int, default=4, help="Adaptive concurrency ceiling (default: 4)")
    args = parser.parse_args()
    if args.max_downloads < 1:
        parser.error("--max-downloads must be positive")
    try:
        accessions = args.accessions
        if args.accession_file:
            accessions += args.accession_file.read_text(encoding="utf-8-sig").split()
        if not accessions:
            parser.error("Provide at least one reference accession")
        accessions = list(dict.fromkeys(map(validate_accession, accessions)))
        directory = args.output.resolve()
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        report = args.report.resolve() if args.report else directory / f"downloads-{stamp}.json"
        directory.mkdir(parents=True, exist_ok=True)
        report.parent.mkdir(parents=True, exist_ok=True)
        with report.open("x", encoding="utf-8"):
            pass
        print(f"Progress report: {report}", flush=True)
        return download_batch(accessions, directory, report, args.max_downloads)
    except (OSError, ValueError) as exc:
        parser.exit(1, f"Download batch stopped: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
