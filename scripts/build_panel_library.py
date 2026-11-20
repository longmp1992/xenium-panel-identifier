#!/usr/bin/env python3
"""Build the bundled reference library of pre-designed 10x Genomics Xenium panels.

Downloads the official ``gene_panel.json`` of every panel listed in
``panel_sources.json`` (links copied from the 10x Genomics Xenium Panel
Designer support pages), parses the gene targets and writes one library file:

    data/xenium_panel_library.json

Nothing is filled in from memory: a panel whose file cannot be downloaded or
parsed is recorded under ``not_retrieved`` with the error and gets no gene list.

Usage:
    python scripts/build_panel_library.py            # build / refresh the library
    python scripts/build_panel_library.py --discover # also look for new panel links

Only the Python standard library is needed.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SOURCES = Path(__file__).resolve().parent / "panel_sources.json"
OUT = ROOT / "data" / "xenium_panel_library.json"
USER_AGENT = "xenium-panel-identifier-library-builder/1.0 (python-urllib)"


def fetch(url: str, retries: int = 3) -> bytes:
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001 - recorded in the library
            last = e
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"{type(last).__name__}: {last}")


def parse_gene_panel(raw: bytes) -> dict:
    """Extract targets from a Xenium gene_panel.json (spec 0.1.0 layout)."""
    doc = json.loads(raw)
    payload = doc.get("payload", doc)
    genes, ensembl, other = [], {}, {}
    for t in payload["targets"]:
        kind = t["type"]["descriptor"]
        data = t["type"].get("data", {})
        if kind == "gene":
            genes.append(data["name"])
            if data.get("id"):
                ensembl[data["name"]] = data["id"]
        else:
            other[kind] = other.get(kind, 0) + 1
    if len(set(genes)) != len(genes):
        raise ValueError("duplicate gene names in panel file")
    panel = payload.get("panel", {})
    ident = panel.get("identity", {})
    return {
        "genes": sorted(genes),
        "ensembl_ids": ensembl,
        "non_gene_target_counts": other,
        "file_panel_name": ident.get("name"),
        "design_id": ident.get("design_id"),
        "design_version": ident.get("version"),
        "file_species": panel.get("species"),
        "file_tissue": panel.get("tissue"),
        "file_num_gene_targets": panel.get("num_gene_targets"),
        "file_chemistry": payload.get("chemistry"),
    }


def discover(known_urls: set[str], index_pages: list[str]) -> None:
    """Best effort: list gene_panel JSON links on the index pages that are not in the manifest."""
    print("\n-- discovery --")
    for page in index_pages:
        try:
            html = fetch(page, retries=1).decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001
            print(f"could not read {page}: {e}\n  (the 10x site rate-limits scripts; open the page in a "
                  "browser and compare its JSON links with panel_sources.json)")
            continue
        found = set(re.findall(r"https://cdn\.10xgenomics\.com/[^\"'\\ ]+?\.json", html))
        new = sorted(found - known_urls)
        gone = sorted(u for u in known_urls if u not in found and "Xenium-panels" in u)
        print(f"{page}\n  json links on page: {len(found)}; not in manifest: {len(new)}")
        for u in new:
            print("  NEW:", u)
        if not found:
            print("  no links found in the HTML (page may be client-rendered) - check it in a browser")
        del gone


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--discover", action="store_true", help="look for panel JSON links missing from the manifest")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()

    src = json.loads(SOURCES.read_text(encoding="utf-8"))
    today = dt.date.today().isoformat()
    panels, failed = [], []
    for p in src["panels"]:
        print(f"{p['id']:<16} ", end="", flush=True)
        try:
            raw = fetch(p["json_url"])
            parsed = parse_gene_panel(raw)
        except Exception as e:  # noqa: BLE001
            print(f"FAILED: {e}")
            failed.append({**p, "error": str(e), "attempted": today})
            continue
        n = len(parsed["genes"])
        warn = []
        if n != p["genes_stated"]:
            warn.append(f"parsed {n} genes but the 10x page states {p['genes_stated']}")
        if parsed["file_num_gene_targets"] not in (None, n):
            warn.append(f"file says num_gene_targets={parsed['file_num_gene_targets']}, parsed {n}")
        if parsed["file_species"] and parsed["file_species"].lower() != p["species"].lower():
            warn.append(f"species mismatch: file={parsed['file_species']}")
        print(f"{n:>5} genes  {parsed['design_id']} v{parsed['design_version']}" + (f"  WARNING: {'; '.join(warn)}" if warn else ""))
        panels.append({
            "id": p["id"],
            "name": p["name"],
            "species": p["species"],
            "chemistry": p["chemistry"],
            "n_genes": n,
            "genes_stated_on_10x_page": p["genes_stated"],
            "source_url": p["json_url"],
            "index_page": p["index_page"],
            "retrieved": today,
            "sha256": hashlib.sha256(raw).hexdigest(),
            "warnings": warn,
            **parsed,
        })

    library = {
        "schema_version": 1,
        "library_version": today,
        "description": "Gene targets of the pre-designed 10x Genomics Xenium panels, parsed from the official gene_panel.json files.",
        "index_pages": src["index_pages"],
        "generated_by": "scripts/build_panel_library.py",
        "panels": panels,
        "not_retrieved": failed,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(library, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"\nwrote {args.out}  ({len(panels)} panels retrieved, {len(failed)} not retrieved)")

    if args.discover:
        discover({p["json_url"] for p in src["panels"]}, src["index_pages"])
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
