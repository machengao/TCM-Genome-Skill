#!/usr/bin/env python3
"""Read the public 2025 Chinese Pharmacopoeia web edition."""

import argparse
from collections import Counter
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import subprocess
import time
from urllib.parse import urlencode


SITE = "https://2025.chp.org.cn/"
BOOK = 1
SECTIONS = {"materia-medica": 9, "formulas": 11}
SCRIPT_VERSION = 1


class ServiceError(Exception):
    def __init__(self, message, status=None):
        super().__init__(message)
        self.status = status


class Paragraphs(HTMLParser):
    def __init__(self, html):
        super().__init__()
        self.rows = []
        self.consuming = False
        self.parts = []
        self.feed(html)
        self.close()

    def handle_starttag(self, tag, attrs):
        if tag == "p":
            self.parts = []
            self.consuming = True
        if tag == "br" and self.consuming:
            self.parts.append(" ")

    def handle_data(self, data):
        if self.consuming:
            self.parts.append(data)

    def handle_endtag(self, tag):
        if tag == "p" and self.consuming:
            self.consuming = False
            text = " ".join("".join(self.parts).split())
            if text:
                self.rows.append(text)


def get_data(resource, **params):
    url = SITE + "front-api/" + resource
    if params:
        url += "?" + urlencode(params)
    for attempt in range(3):
        completed = subprocess.run(
            ["curl", "--silent", "--show-error", "--max-time", "15", "--write-out", "\n%{http_code}", url],
            capture_output=True, check=False,
        )
        body, _, http_status = completed.stdout.rpartition(b"\n")
        status = int(http_status) if http_status.isdigit() else 0
        if completed.returncode == 0 and status == 200:
            result = json.loads(body)
            if result.get("code") != 200 or result.get("data") is None:
                code = result.get("code")
                raise ServiceError(f"Unexpected service response for {resource}: {code}", code)
            return result["data"]
        if status not in (0, 408, 500, 502, 503, 504) or attempt == 2:
            detail = completed.stderr.decode("utf-8", errors="replace").strip()
            raise ServiceError(f"{resource}: HTTP {status or 'unavailable'} {detail}".strip(), status)
        time.sleep(attempt + 1)


def list_entries(section):
    directory_id = SECTIONS[section]
    result = get_data("children", bookId=BOOK, directoryId=directory_id)
    entries = [row for row in result["children"] if row["type"] == "entry"]
    if result["directory"]["directoryId"] != directory_id or not entries:
        raise ValueError("The official directory failed validation")
    if len({row["id"] for row in entries}) != len(entries):
        raise ValueError("The official directory has duplicate entry IDs")
    return entries


def read_entry(stub, section):
    entry_id = stub["id"]
    row = get_data(f"entry/{entry_id}")
    if (row["bookId"], row["directoryId"], row["entryId"], row["title"]) != (
        BOOK, SECTIONS[section], entry_id, stub["title"]
    ):
        raise ValueError("The official entry does not match the directory")
    paragraphs = Paragraphs(row["htmlContent"]).rows
    if not paragraphs:
        raise ValueError("No readable paragraphs in the official entry")

    heading = re.compile(r"^【[^】]+】")
    if section == "formulas":
        start = next((i for i, line in enumerate(paragraphs) if line.startswith("【处方】")), None)
        lines = [] if start is None else [paragraphs[start].removeprefix("【处方】").strip()]
        following = [] if start is None else paragraphs[start + 1:]
    else:
        lines = []
        following = paragraphs
    for line in following:
        if heading.match(line):
            break
        lines.append(line)
    return {
        "status": "read",
        "title": row["title"],
        "webPage": row["pageNum"],
        "source": f"{SITE}#/item?bookId={BOOK}&entryId={entry_id}",
        "sectionText": [line for line in lines if line],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("section", choices=SECTIONS)
    parser.add_argument("--output", type=Path, required=True, help="JSON report; reused to resume incomplete runs")
    parser.add_argument("--title", action="append", help="Exact entry title; omit to read the whole section")
    parser.add_argument("--interval", type=float, default=0.35, help="Minimum seconds between entry requests")
    args = parser.parse_args()
    try:
        catalog = list_entries(args.section)
    except (OSError, ServiceError, ValueError, KeyError, TypeError) as exc:
        parser.exit(1, f"Cannot read official directory: {exc}\n")
    selected = [row for row in catalog if args.title is None or row["title"] in args.title]
    missing_titles = sorted(set(args.title or ()) - {row["title"] for row in selected})
    report = {"schema": SCRIPT_VERSION, "edition": 2025, "book": BOOK, "section": args.section, "entries": {}}
    if args.output.exists():
        previous = json.loads(args.output.read_text(encoding="utf-8"))
        if any(previous.get(key) != report[key] for key in ("schema", "edition", "book", "section")):
            parser.error("The output file belongs to another edition, section or schema")
        report["entries"] = previous["entries"]
    report["catalog"] = catalog
    report["selection"] = [row["id"] for row in selected]
    report["missingTitles"] = missing_titles
    args.output.parent.mkdir(parents=True, exist_ok=True)

    def checkpoint():
        temporary = args.output.with_name(args.output.name + ".partial")
        temporary.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
        temporary.replace(args.output)

    try:
        for number, stub in enumerate(selected, 1):
            key = str(stub["id"])
            cached = report["entries"].get(key, {})
            if cached.get("status") == "read" and cached.get("title") == stub["title"]:
                continue
            try:
                report["entries"][key] = read_entry(stub, args.section)
            except (ServiceError, ValueError, KeyError, TypeError) as exc:
                report["entries"][key] = {"status": "error", "title": stub["title"], "error": str(exc)}
                if isinstance(exc, ServiceError) and exc.status in (401, 403, 429):
                    print(f"Stopped after service status {exc.status}.")
                    break
            if number % 50 == 0:
                checkpoint()
                print(f"{number}/{len(selected)} entries processed", flush=True)
            time.sleep(max(0.25, args.interval))
    finally:
        checkpoint()

    statuses = Counter(report["entries"].get(str(row["id"]), {}).get("status", "pending") for row in selected)
    print(json.dumps({"catalog": len(catalog), "selected": len(selected), **statuses, "missingTitles": missing_titles},
                     ensure_ascii=False))
    return int(bool(missing_titles) or statuses["read"] != len(selected))


if __name__ == "__main__":
    raise SystemExit(main())
