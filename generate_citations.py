#!/usr/bin/env python3
"""
Fetch publications from ORCID + CrossRef and write citations.yaml.

Outputs:
  - citations.yaml          (webapp root, read by app.py)
  - static/lab_webpage/citations.yaml  (static copy)

Run:  python3 generate_citations.py
"""

import json
import time
import yaml
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import URLError

ORCID = "0000-0002-1145-6199"

HERE = Path(__file__).resolve().parent
OUT_PATHS = [
    HERE / "citations.yaml",
    HERE / "static" / "lab_webpage" / "citations.yaml",
]

CROSSREF_HEADERS = {
    "User-Agent": "citations-generator/1.0 (mailto:den.torre.94@gmail.com)"
}


def fetch_json(url, headers=None):
    req = Request(url, headers=headers or {})
    with urlopen(req, timeout=15) as resp:
        return json.loads(resp.read())


def get_orcid_dois(orcid):
    url = f"https://pub.orcid.org/v3.0/{orcid}/works"
    data = fetch_json(url, {"Accept": "application/json"})
    dois = []
    for work in data.get("group", []):
        ids = []
        for summary in work.get("work-summary", []):
            ids.extend(summary.get("external-ids", {}).get("external-id", []))
        for ext_id in ids:
            if (
                ext_id.get("external-id-type") == "doi"
                and ext_id.get("external-id-relationship")
                in ("self", "version-of", "part-of")
            ):
                dois.append(ext_id["external-id-value"])
                break
    # deduplicate while preserving order
    seen = set()
    return [d for d in dois if not (d in seen or seen.add(d))]


def crossref_citation(doi, orcid):
    url = f"https://api.crossref.org/works/{doi}"
    try:
        data = fetch_json(url, CROSSREF_HEADERS)
    except (URLError, Exception) as e:
        print(f"    CrossRef error: {e}")
        return None

    msg = data.get("message", {})

    title = (msg.get("title") or [""])[0].strip()

    authors = []
    for author in msg.get("author", []):
        given = author.get("given", "").strip()
        family = author.get("family", "").strip()
        name = f"{given} {family}".strip()
        if name:
            authors.append(name)

    container = (msg.get("container-title") or [""])[0].strip()
    publisher = container or msg.get("publisher", "").strip()

    _preprint_publishers = {"10.1101/": "bioRxiv", "10.64898/": "openRxiv"}
    for prefix, name in _preprint_publishers.items():
        if doi.startswith(prefix):
            publisher = name
            break

    date_str = ""
    for date_field in ("published", "published-print", "published-online", "issued"):
        parts = (msg.get(date_field) or {}).get("date-parts", [[]])
        if parts and parts[0]:
            p = parts[0]
            year = p[0] if len(p) > 0 else None
            month = p[1] if len(p) > 1 else 1
            day = p[2] if len(p) > 2 else 1
            if year:
                date_str = f"{year}-{str(month).zfill(2)}-{str(day).zfill(2)}"
                break

    link = msg.get("URL", "").strip() or f"https://doi.org/{doi}"

    return {
        "id": f"doi:{doi}",
        "title": title,
        "authors": authors,
        "publisher": publisher,
        "date": date_str,
        "link": link,
        "orcid": orcid,
    }


def main():
    print(f"Fetching works for ORCID {ORCID}…")
    dois = get_orcid_dois(ORCID)
    print(f"Found {len(dois)} DOIs\n")

    citations = []
    for i, doi in enumerate(dois, 1):
        print(f"  [{i:>2}/{len(dois)}] {doi}")
        citation = crossref_citation(doi, ORCID)
        if citation and citation.get("date"):
            citations.append(citation)
        time.sleep(0.1)

    citations.sort(key=lambda c: c["date"], reverse=True)
    print(f"\n{len(citations)} citations generated")

    for path in OUT_PATHS:
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("# DO NOT EDIT, GENERATED AUTOMATICALLY\n\n")
            yaml.dump(
                citations, fh,
                default_flow_style=False,
                sort_keys=False,
                allow_unicode=True,
            )
        print(f"Written → {path}")


if __name__ == "__main__":
    main()
