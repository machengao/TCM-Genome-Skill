#!/usr/bin/env python3
"""Search public genome metadata for a list of checked pharmacopoeial taxa."""

import argparse
from collections import Counter
from datetime import datetime, timezone
from html import unescape
import json
from pathlib import Path
import re
import subprocess
import time
from urllib.parse import quote, urlencode


NCBI = "https://api.ncbi.nlm.nih.gov/datasets/v2"
GWH = "https://ngdc.cncb.ac.cn/gwh"
ENA = "https://www.ebi.ac.uk/ena/portal/api/search"
JGI = "https://files.jgi.doe.gov"
PHYCOCOSM = "https://phycocosm.jgi.doe.gov"
GWH_COLUMNS = (
    "", "scientificName", "commonNames", "group", "source", "accession",
    "refseqAccession", "hasReannotation", "genomeRepresentation", "assemblyLevel",
    "genomeSize", "chrCount", "gcContent", "releaseDate", "sourceReleaseDate",
    "dna", "gff", "rna", "cds", "protein", "refseqDna", "refseqGff", "refseqRna",
    "refseqCds", "refseqProtein",
)
JGI_SOURCES = {
    "mycocosm": ("mycocosm_file_list", "assembly"),
    "phytozome": ("phytozome_file_list", "fasta"),
    "phycocosm": ("phycocosm_file_list", ""),
    "img": ("img_file_list", ""),
}
BASELINE = ("ncbi", "gwh")
SUPPLEMENTS = {
    "fungi": ("mycocosm", "ena"),
    "plants": ("phytozome", "ena"),
    "animals": ("ena",),
    "algae": ("phycocosm", "ena"),
    "microbes": ("img", "ena"),
    "other": ("ena",),
}
ALIASES = {
    "ncbi datasets": "ncbi", "genbank": "ncbi", "refseq": "ncbi",
    "ngdc": "gwh", "cncb": "gwh", "genome warehouse": "gwh",
    "european nucleotide archive": "ena", "jgi mycocosm": "mycocosm",
    "jgi phytozome": "phytozome", "jgi phycocosm": "phycocosm", "jgi img": "img",
}
ACCESSION = re.compile(r"(?:GC[AF]_\d+|GWH[A-Z]+\d+)(?:\.\d+)?\Z")
SEQUENCE_TYPES = ("nuclear", "mitochondrion", "chloroplast")
SEQUENCE_ACCESSION = re.compile(r"[A-Z][A-Z_]*\d+\.\d+\Z")


class QueryError(Exception):
    def __init__(self, detail, status="error"):
        super().__init__(detail)
        self.status = status


def request(url, params=None, payload=None, accept="application/json"):
    if params:
        url += "?" + urlencode(params)
    args = [
        "curl", "--silent", "--show-error", "--connect-timeout", "8", "--max-time",
        "25", "--max-filesize", "12000000", "--header", f"Accept: {accept}",
        "--write-out", "\n%{http_code}", url,
    ]
    body = None
    if payload is not None:
        args.extend(["--header", "Content-Type: application/json", "--data-binary", "@-"])
        body = json.dumps(payload).encode()
    for attempt in range(2):
        result = subprocess.run(args, input=body, capture_output=True, check=False)
        data, _, footer = result.stdout.rpartition(b"\n")
        code = int(footer) if footer.isdigit() else 0
        if result.returncode == 0 and code == 200:
            return data
        if code in (401, 403):
            raise QueryError(f"HTTP {code}", "access_required")
        if code == 429:
            raise QueryError("HTTP 429", "rate_limited")
        if attempt or code not in (0, 408, 500, 502, 503, 504) or result.returncode == 63:
            detail = result.stderr.decode("utf-8", errors="replace").strip()
            raise QueryError(f"HTTP {code or 'unavailable'} {detail}".strip())
        time.sleep(1)


def get_json(url, params=None, payload=None):
    data = request(url, params, payload)
    try:
        return json.loads(data)
    except ValueError as exc:
        raise QueryError("The service did not return JSON") from exc


def accession(value):
    return value if isinstance(value, str) and ACCESSION.fullmatch(value) else None


def record(source, identifier, organism, url, **metadata):
    return dict(source=source, accession=accession(identifier), organism=organism, url=url, **metadata)


def match_name(organism, names):
    if not organism:
        return False
    label = " ".join(organism.lower().split())
    words = [" ".join(term.lower().split()) for term in names]
    return any(label == term or label.startswith(term + " ") for term in words)


def ncbi_search(term, names):
    url = f"{NCBI}/genome/taxon/{quote(term, safe='')}/dataset_report"
    results = []
    token = None
    for _ in range(100):
        parameters = {"page_size": 100}
        if token:
            parameters["page_token"] = token
        data = get_json(url, parameters)
        if not isinstance(data, dict) or data.get("message") or (
            "reports" not in data and data and data.get("total_count") != 0
        ):
            raise QueryError("Unexpected NCBI dataset report")
        for item in data.get("reports", []):
            identifier = accession(item["accession"])
            name = item["organism"]["organism_name"]
            if not identifier:
                raise QueryError("NCBI returned an invalid assembly accession")
            info = item.get("assembly_info", {})
            organism = item["organism"]
            results.append(record(
                "ncbi", identifier, name, f"https://www.ncbi.nlm.nih.gov/datasets/genome/{identifier}/",
                assembly_level=info.get("assembly_level"), assembly_name=info.get("assembly_name"),
                tax_id=organism.get("tax_id"), length_bp=item.get("assembly_stats", {}).get("total_sequence_length"),
                strain=organism.get("infraspecific_names", {}).get("strain"),
                paired_accession=accession(item.get("paired_accession")),
                biosample=info.get("biosample", {}).get("accession"), release_date=info.get("release_date"),
                sequence_type="assembly_unverified",
                match="matched" if match_name(name, names) else "candidate",
            ))
        next_token = data.get("next_page_token")
        if not next_token:
            return results
        if token == next_token:
            break
        token = next_token
    raise QueryError("NCBI pagination did not complete")


def ncbi_organelle_search(term, names, wanted):
    url = f"{NCBI}/organelle/taxon/{quote(term, safe='')}/dataset_report"
    results = []
    token = None
    for _ in range(100):
        parameters = {"page_size": 100}
        if token:
            parameters["page_token"] = token
        data = get_json(url, parameters)
        if not isinstance(data, dict) or (
            "reports" not in data and data and data.get("total_count") != 0
        ):
            raise QueryError("Unexpected NCBI organelle report")
        for item in data.get("reports", []):
            description = item["description"]
            sequence_type = description.strip().lower()
            if sequence_type not in wanted:
                continue
            ids = {key: (item.get(key) or {}).get("accession") for key in ("refseq", "genbank")}
            ids = {key: value for key, value in ids.items()
                   if isinstance(value, str) and SEQUENCE_ACCESSION.fullmatch(value)}
            identifier = ids.get("refseq") or ids.get("genbank")
            if not identifier:
                raise QueryError("NCBI organelle record has no valid sequence accession")
            organism = item["organism"]
            name = organism["organism_name"]
            results.append(record(
                "ncbi", None, name, f"https://www.ncbi.nlm.nih.gov/nuccore/{identifier}",
                record_id=identifier, sequence_accessions=ids, sequence_type=sequence_type,
                type_evidence=description, length_bp=item.get("length"), topology=item.get("topology"),
                tax_id=organism.get("tax_id"), sample_names=organism.get("infraspecific_names"),
                match="matched" if match_name(name, names) else "candidate",
            ))
        next_token = data.get("next_page_token")
        if not next_token:
            return results
        if next_token == token:
            break
        token = next_token
    raise QueryError("NCBI organelle pagination did not complete")


def gwh_sequence_type(url, representation):
    label = representation or ""
    if label.strip().lower() == "full genome":
        return "assembly_unverified", representation
    if label.strip().lower() == "partial genome":
        html = request(url, accept="text/html").decode("utf-8")
        field = re.search(r"<dt>\s*Genome Representation\s*</dt>\s*<dd>(.*?)</dd>", html, re.S | re.I)
        if not field:
            raise QueryError("Missing GWH genome representation detail")
        label = " ".join(unescape(re.sub(r"<[^>]*>", " ", field[1])).split())
    kinds = re.findall(r"\b(chloroplast|mitochondri(?:on|al|a))\b", label.lower())
    kinds = {"chloroplast" if name == "chloroplast" else "mitochondrion" for name in kinds}
    return (next(iter(kinds)) if len(kinds) == 1 else "unverified"), label


def gwh_search(term, names, wanted=SEQUENCE_TYPES):
    columns = [
        {"data": column, "name": column, "searchable": bool(column), "orderable": bool(column),
         "search": {"value": term if column == "scientificName" else "", "regex": False}}
        for column in GWH_COLUMNS
    ]
    results = []
    start = 0
    while True:
        data = get_json(
            GWH + "/api/esdt/gwh_assemblies/query",
            payload={"draw": 1, "start": start, "length": 100, "columns": columns,
                     "order": [{"column": 5, "dir": "asc"}], "search": {"value": "", "regex": False}},
        )
        rows = data["data"]
        total = int(data["recordsFiltered"])
        if not isinstance(rows, list) or total > 10000 or start + len(rows) > total:
            raise QueryError("Unexpected or overly broad GWH directory result")
        for item in rows:
            identifier = accession(item["accession"])
            name = item["scientificName"]
            if not identifier:
                raise QueryError("GWH returned an invalid assembly accession")
            matched = match_name(name, names)
            if not matched and not match_name(name, [term]):
                continue
            direct = item.get("source") == "Direct submission"
            target = f"Assembly/{item['primaryId']}/show" if direct else f"Genome/{item['genomeId']}/show"
            url = f"{GWH}/{target}"
            details = {}
            representation = item.get("genomeRepresentation")
            try:
                if not direct and representation and representation.strip().lower() == "partial genome":
                    # The mirror links to a species page, which cannot establish an assembly's molecule type.
                    if "nuclear" not in wanted:
                        continue
                    sequence_type, evidence = "unverified", representation
                else:
                    sequence_type, evidence = gwh_sequence_type(url, representation)
            except QueryError as exc:
                sequence_type, evidence = "unverified", representation
                details["type_error"] = str(exc)
            if sequence_type in SEQUENCE_TYPES[1:] and sequence_type not in wanted:
                continue
            if "nuclear" not in wanted and sequence_type not in wanted and sequence_type != "unverified":
                continue
            if direct and identifier.startswith("GWH"):
                detail_url = GWH + "/api/public/assembly/" + identifier
                try:
                    detail = get_json(detail_url)
                    if detail.get("assemblyAccession") != identifier:
                        raise QueryError("GWH detail did not verify the indexed accession")
                    details.update(tax_id=detail.get("taxId"), biosample=detail.get("biosampleAccession"))
                except (QueryError, AttributeError) as exc:
                    details["detail_error"] = str(exc)
                details["detail_url"] = detail_url
            results.append(record(
                "gwh", identifier, name, url, sequence_type=sequence_type, type_evidence=evidence,
                gwh_source=item.get("source"), assembly_name=item.get("assemblyName"),
                assembly_level=item.get("assemblyLevel"), release_date=item.get("sourceReleaseDate"),
                length_mb=item.get("genomeSize"), paired_accession=accession(item.get("refseqAccession")),
                match="matched" if matched and sequence_type != "unverified" else "candidate", **details,
            ))
        start += len(rows)
        if start == total:
            return results
        if not rows:
            raise QueryError("GWH pagination ended before the reported total")


def ena_search(term, names):
    # NCBI taxonomy IDs include assemblies registered below species rank.
    expression = f"tax_tree({term})" if term.isascii() and term.isdecimal() else f'scientific_name="{term}"'
    data = get_json(ENA, {"result": "assembly", "query": expression, "format": "json", "limit": 0,
                         "fields": "accession,version,scientific_name,tax_id,assembly_level,strain,sample_accession"})
    if not isinstance(data, list):
        raise QueryError("Unexpected ENA assembly result")
    results = []
    for item in data:
        identifier = item["accession"]
        version = item.get("version")
        if accession(identifier) and "." not in identifier and str(version).isdigit():
            identifier += f".{version}"
        name = item["scientific_name"]
        results.append(record(
            "ena", identifier, name, "https://www.ebi.ac.uk/ena/browser/view/" + quote(identifier, safe=""),
            tax_id=item.get("tax_id"), assembly_level=item.get("assembly_level"),
            strain=item.get("strain"), biosample=item.get("sample_accession"),
            record_id=identifier, match="matched" if match_name(name, names) else "candidate",
        ))
    return results


def phycocosm_catalog(term, names):
    html = request(PHYCOCOSM + "/phycocosm/home/releases", {"flt": term, "pub": "public"}, accept="text/html")
    table = re.search(r'<table\b[^>]*\bid="catalog"[^>]*>(.*?)</table>', html.decode("utf-8"), re.S)
    if not table:
        raise QueryError("PhycoCosm did not return its published genome catalog")
    rows = re.findall(r'<td\b[^>]*><a href="(/[\w.-]+)">([^<]+)</a></td>', table[1])
    return [
        record("phycocosm", None, unescape(label), PHYCOCOSM + path, dataset_id=path[1:],
               match="candidate", note="Public PhycoCosm catalog; verify cross-references before merging.")
        for path, label in rows if match_name(unescape(label), [*names, term])
    ]


def jgi_search(source, term, names):
    endpoint, file_type = JGI_SOURCES[source]
    url = f"{JGI}/{endpoint}/"
    results = []
    page = 1
    while page <= 20:
        params = {"api_version": 2, "q": f'"{term}"', "p": page, "x": 3}
        if file_type:
            params["ff[file_type]"] = file_type
        data = get_json(url, params)
        organisms = data["organisms"]
        if not isinstance(organisms, list) or "next_page" not in data:
            raise QueryError("Unexpected JGI metadata result")
        for dataset in organisms:
            if dataset["kingdom"] != source or dataset["label"] != "Genome":
                continue
            metadata = dataset["top_hit"]["metadata"]
            labels = (
                metadata.get("ncbi_taxon", {}).get("ncbi_taxon_species"),
                metadata.get("phytozome", {}).get("proteome_name"),
                metadata.get("gold_data", {}).get("organism_name"),
                metadata.get("portal_metadata", {}).get("long_display_name"),
                metadata.get("analysis_project", {}).get("analysis_project_name"),
            )
            name = next((label for label in labels if match_name(label, [*names, term])), None)
            if not name:
                continue
            dataset_id = str(dataset["id"])
            results.append(record(
                source, None, name, "https://data.jgi.doe.gov/search?" + urlencode({"q": term}),
                dataset_id=dataset_id, file_name=dataset["top_hit"].get("file_name"),
                portal_id=metadata.get(f"{source}_portal_id"), match="candidate",
                note="JGI genome dataset. Verify assembly version and cross-references before merging.",
            ))
        next_page = data["next_page"]
        if not next_page:
            if not results and source == "phycocosm":
                return phycocosm_catalog(term, names)
            return results
        if type(next_page) is not int or next_page <= page:
            break
        page = next_page
    raise QueryError("JGI pagination did not complete")


def search_source(source, terms, names, attempted, wanted=("nuclear",)):
    has_nuclear = "nuclear" in wanted
    organelles = set(wanted) - {"nuclear"}
    if source == "ncbi" and has_nuclear and organelles:
        assemblies, assembly_log = search_source(source, terms, names, attempted, ("nuclear",))
        sequences, sequence_log = search_source(source, terms, names, attempted, tuple(sorted(organelles)))
        return assemblies + sequences, assembly_log + sequence_log
    if source not in (*BASELINE, "ena", *JGI_SOURCES):
        return [], [{"status": "manual_required", "detail": "No automated adapter for this source"}]
    if not has_nuclear and source not in BASELINE:
        return [], [{"status": "manual_required", "detail": "No organelle adapter for this source"}]
    results, audit = [], []
    for term in terms:
        try:
            key = (source, term, tuple(names), tuple(wanted))
            if key not in attempted:
                if source == "ncbi":
                    rows = ncbi_search(term, names) if has_nuclear else []
                    if organelles:
                        rows.extend(ncbi_organelle_search(term, names, organelles))
                elif source == "gwh":
                    rows = gwh_search(term, names, wanted)
                elif source == "ena":
                    rows = ena_search(term, names)
                else:
                    rows = jgi_search(source, term, names)
                attempted[key] = rows
                time.sleep(0.25)
            rows = attempted[key]
            matched = sum(row["match"] == "matched" for row in rows)
            status = "hit" if matched else "candidate" if rows else "no_hit"
            audit.append({"term": term, "status": status, "matched": matched, "candidates": len(rows) - matched,
                          "types_found": dict(Counter(row.get("sequence_type", "assembly_unverified") for row in rows))})
            results.extend(rows)
        except (OSError, QueryError, KeyError, TypeError, ValueError, AttributeError) as exc:
            status = exc.status if isinstance(exc, QueryError) else "error"
            audit.append({"term": term, "status": status, "detail": str(exc)})
            if status in ("access_required", "rate_limited"):
                audit.extend({"term": pending, "status": "not_run"} for pending in terms[len(audit):])
                break
    if organelles and source not in BASELINE:
        audit.append({"status": "manual_required", "detail": "No organelle adapter for this source",
                      "sequence_types": sorted(organelles)})
    return results, [dict({"sequence_types": list(wanted)}, **result) for result in audit]


def combine_records(records):
    groups = []
    matches = {}
    for row in records:
        keys = {value for value in (row["accession"], row.get("paired_accession"))
                if value and (value.startswith("GWH") or "." in value)}
        if not keys:
            keys.add((row["source"], row["accession"] or row.get("dataset_id") or row.get("record_id") or row["url"]))
        existing = [matches[key] for key in keys if key in matches]
        group = existing[0] if existing else {"keys": set(), "evidence": []}
        if not existing:
            groups.append(group)
        for other in existing[1:]:
            if other is not group:
                group["keys"].update(other["keys"])
                group["evidence"].extend(other["evidence"])
                other["keys"].clear()
        group["keys"].update(keys)
        for key in group["keys"]:
            matches[key] = group
        if row not in group["evidence"]:
            group["evidence"].append(row)
    return [
        {"accessions": sorted({value for row in group["evidence"]
                               for value in (row["accession"], row.get("paired_accession")) if value}),
         "sources": sorted({row["source"] for row in group["evidence"]}),
         "sequence_types": sorted({row.get("sequence_type", "assembly_unverified") for row in group["evidence"]}),
         "status": "matched" if any(row["match"] == "matched" for row in group["evidence"]) else "candidate",
         "evidence": group["evidence"]}
        for group in groups if group["keys"]
    ]


def source_name(value):
    name = " ".join(value.lower().split())
    return ALIASES.get(name, name)


def lookup(taxon, force_sources, attempted, wanted=("nuclear",)):
    names = list(dict.fromkeys([taxon["name"], *taxon.get("aliases", [])]))
    candidates = list(dict.fromkeys(taxon.get("candidate_names", [])))
    preferred = list(dict.fromkeys(map(source_name, [*force_sources, *taxon.get("required_sources", [])])))
    required = list(dict.fromkeys([*BASELINE, *preferred]))
    records, sources = [], {}

    def check(source, mandatory):
        terms = [*names, *candidates]
        if taxon.get("tax_id") and source in ("ncbi", "ena"):
            terms.insert(0, str(taxon["tax_id"]))
        rows, log = search_source(source, terms, names, attempted, wanted)
        records.extend(rows)
        sources[source] = {"required": mandatory, "queries": log}

    for source in required:
        check(source, True)
    if "nuclear" in wanted and not any(
        row["match"] == "matched" and row.get("sequence_type") not in SEQUENCE_TYPES[1:]
        for row in records if row["source"] in BASELINE
    ):
        for source in SUPPLEMENTS[taxon["group"]]:
            if source not in sources:
                check(source, False)
    completed = {"hit", "no_hit", "candidate"}
    incomplete = any(q["status"] not in completed for report in sources.values() for q in report["queries"])
    groups = combine_records(records)
    status = (
        "hit" if any(group["status"] == "matched" for group in groups)
        else "candidate" if groups else "incomplete" if incomplete else "no_hit"
    )
    matched_types = {row.get("sequence_type", "assembly_unverified") for row in records
                     if row["match"] == "matched"}
    if "assembly_unverified" in matched_types:
        matched_types.add("nuclear")
    return {"taxon": taxon, "status": status, "sources": sources, "assemblies": groups,
            "literature_followup": [kind for kind in wanted if kind not in matched_types],
            "required_complete": all(
                q["status"] in completed for report in sources.values() if report["required"] for q in report["queries"]
            ), "coverage_complete": not incomplete}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--require", action="append", default=[], help="Additional required source (repeatable)")
    parser.add_argument("--sequence-type", action="append", choices=SEQUENCE_TYPES,
                        help="Requested sequence type (repeatable); defaults to nuclear")
    args = parser.parse_args()
    taxa = json.loads(args.input.read_text(encoding="utf-8-sig"))
    if not isinstance(taxa, list) or not taxa:
        parser.error("Expected a nonempty list of taxa")
    for item in taxa:
        if not isinstance(item, dict) or item.get("group") not in SUPPLEMENTS:
            parser.error(f"Each taxon requires a group from: {', '.join(SUPPLEMENTS)}")
        if not isinstance(item.get("name"), str) or not item["name"].strip():
            parser.error("Each taxon requires a scientific name")
        for key in ("aliases", "candidate_names", "required_sources"):
            values = item.get(key, [])
            if not isinstance(values, list) or any(not isinstance(value, str) or not value.strip() for value in values):
                parser.error(f"{key} must be a list of nonempty strings")
        if any('"' in term or "\\" in term for term in [item["name"], *item.get("aliases", []), *item.get("candidate_names", [])]):
            parser.error("Scientific names cannot contain ASCII double quotes or backslashes")
        if item.get("tax_id") and not re.fullmatch(r"[1-9][0-9]*", str(item["tax_id"])):
            parser.error("tax_id must be a positive NCBI taxonomy ID")
    if args.input.resolve() == args.output.resolve():
        parser.error("Input and output must be different files")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    attempted, results = {}, []
    wanted = tuple(dict.fromkeys(args.sequence_type or ["nuclear"]))
    for item in taxa:
        result = lookup(item, args.require, attempted, wanted)
        results.append(result)
        print(f"{item['name']}: {result['status']}", flush=True)
    report = {"queried_at": datetime.now(timezone.utc).isoformat(), "sequence_types": wanted, "results": results}
    temp = args.output.with_name(args.output.name + ".partial")
    temp.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    temp.replace(args.output)
    print(json.dumps(dict(Counter(row["status"] for row in results))))
    return int(any(not row["required_complete"] for row in results))


if __name__ == "__main__":
    raise SystemExit(main())
