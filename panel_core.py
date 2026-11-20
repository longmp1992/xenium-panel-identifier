"""Core logic: read the feature list of a Xenium output file and match it to a panel.

No web framework imports here, so it can be used from the app, the command line
helper and tests alike.
"""
from __future__ import annotations

import csv
import gzip
import io
import json
from dataclasses import dataclass, field
from pathlib import Path

LIBRARY_PATH = Path(__file__).resolve().parent / "data" / "xenium_panel_library.json"

# Name prefixes of non-gene features written by Xenium Onboard Analysis.
CONTROL_PREFIXES = {
    "NegControlProbe_": "negative_control_probe",
    "NegControlCodeword_": "negative_control_codeword",
    "UnassignedCodeword_": "unassigned_codeword",
    "DeprecatedCodeword_": "deprecated_codeword",
    "BLANK_": "blank_codeword",
    "Intergenic_Region_": "genomic_control_probe",
    "antisense_": "negative_control_probe",
}

# Calls
MATCH_MIN = 0.95      # >= this fraction of a panel's genes present -> the panel was used
PARTIAL_MIN = 0.50    # between PARTIAL_MIN and MATCH_MIN -> reported as uncertain


@dataclass
class FeatureTable:
    """Distinct features of one dataset."""
    source_kind: str
    # feature name -> {"category": str|None, "is_gene": bool|None, "count": int|None}
    features: dict = field(default_factory=dict)
    n_rows: int | None = None
    notes: list = field(default_factory=list)


# --------------------------------------------------------------------------- readers

def read_parquet_features(path, progress=None) -> FeatureTable:
    """Distinct feature_name values (+ category columns, + transcript counts).

    Only the needed columns are read, one row group at a time, so memory stays
    small even for multi-GB transcripts.parquet files.
    """
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq

    with open(path, "rb") as fh:  # file handle: also works for UNC / network paths
        pf = pq.ParquetFile(fh)
        names = set(pf.schema_arrow.names)
        if "feature_name" not in names:
            raise ValueError(
                "This parquet file has no 'feature_name' column, so it does not look like a Xenium "
                f"transcripts.parquet. Columns found: {', '.join(pf.schema_arrow.names)}"
            )
        cols = ["feature_name"] + [c for c in ("codeword_category", "is_gene") if c in names]
        ft = FeatureTable("transcripts.parquet", n_rows=pf.metadata.num_rows)
        agg: dict = {}
        n_rg = pf.metadata.num_row_groups
        for i in range(n_rg):
            tbl = pf.read_row_group(i, columns=cols)
            casts = []
            for c in cols:
                col = tbl[c]
                if pa.types.is_dictionary(col.type):
                    col = pc.cast(col, col.type.value_type)
                if pa.types.is_binary(col.type) or pa.types.is_large_binary(col.type):
                    col = pc.cast(col, pa.string())  # older XOA wrote bytes
                casts.append(col)
            tbl = pa.table(casts, names=cols)
            grouped = tbl.group_by(cols).aggregate([([], "count_all")]).to_pydict()
            for j, name in enumerate(grouped["feature_name"]):
                if name is None:
                    continue
                rec = agg.setdefault(name, {"category": None, "is_gene": None, "count": 0})
                rec["count"] += grouped["count_all"][j]
                if "codeword_category" in grouped and grouped["codeword_category"][j] is not None:
                    rec["category"] = grouped["codeword_category"][j]
                if "is_gene" in grouped and grouped["is_gene"][j] is not None:
                    rec["is_gene"] = bool(grouped["is_gene"][j])
            if progress:
                progress((i + 1) / n_rg)
        ft.features = agg
        ft.notes.append(f"columns used: {', '.join(cols)}; writer: {pf.metadata.created_by}")
    return ft


def _open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8-sig", newline="")
    return open(path, encoding="utf-8-sig", newline="")


def read_gene_panel_json(path) -> FeatureTable:
    """A Xenium output ``gene_panel.json``."""
    with open(path, encoding="utf-8") as fh:
        doc = json.load(fh)
    payload = doc.get("payload", doc)
    if "targets" not in payload:
        raise ValueError("JSON file has no 'targets' list - not a Xenium gene_panel.json")
    ft = FeatureTable("gene_panel.json")
    for t in payload["targets"]:
        kind = t["type"]["descriptor"]
        name = t["type"].get("data", {}).get("name")
        if name:
            ft.features[name] = {"category": kind, "is_gene": kind == "gene", "count": None}
    ident = payload.get("panel", {}).get("identity", {})
    if ident:
        ft.notes.append(f"panel identity recorded in the file: {ident.get('name')} "
                        f"(design_id {ident.get('design_id')}, version {ident.get('version')})")
    return ft


def read_feature_list(path) -> FeatureTable:
    """Tabular / plain lists: output of extract_features.py, 10x features.tsv(.gz), one name per line."""
    with _open_text(path) as fh:
        text = fh.read()
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    if not lines:
        raise ValueError("The file is empty.")
    delim = "\t" if "\t" in lines[0] else ("," if "," in lines[0] else None)
    rows = list(csv.reader(lines, delimiter=delim)) if delim else [[ln.strip()] for ln in lines]
    header = [h.strip().lower() for h in rows[0]]
    ft = FeatureTable("feature list")

    if "feature_name" in header:  # extract_features.py output
        ix = {h: k for k, h in enumerate(header)}
        for r in rows[1:]:
            def get(col):
                k = ix.get(col)
                return r[k].strip() if k is not None and k < len(r) and r[k].strip() != "" else None
            isg = get("is_gene")
            cnt = get("n_transcripts")
            ft.features[r[ix["feature_name"]].strip()] = {
                "category": get("codeword_category"),
                "is_gene": None if isg is None else isg.lower() in ("true", "1"),
                "count": int(cnt) if cnt and cnt.isdigit() else None,
            }
    elif len(rows[0]) >= 3 and rows[0][0].upper().startswith(("ENS", "NEGCONTROL", "BLANK", "UNASSIGNED")):
        # 10x features.tsv: id, name, feature type (no header)
        ft.source_kind = "features.tsv (cell_feature_matrix)"
        for r in rows:
            ftype = r[2].strip()
            ft.features[r[1].strip()] = {"category": ftype, "is_gene": ftype == "Gene Expression", "count": None}
    else:  # first column = names; drop an obvious header
        start = 1 if header[0] in ("gene", "genes", "feature", "features", "name", "gene_name", "symbol") else 0
        for r in rows[start:]:
            if r and r[0].strip():
                ft.features[r[0].strip()] = {"category": None, "is_gene": None, "count": None}
    return ft


def read_any(path, progress=None) -> FeatureTable:
    """Dispatch on file content / extension."""
    p = str(path)
    with open(p, "rb") as fh:
        magic = fh.read(4)
    if magic == b"PAR1":
        return read_parquet_features(p, progress)
    low = p.lower()
    if low.endswith((".parquet", ".pq")):
        raise ValueError("The file has a .parquet name but is not a valid Parquet file (truncated download?).")
    if low.endswith((".json", ".xenium")):
        return read_gene_panel_json(p)
    return read_feature_list(p)


# --------------------------------------------------------------------------- classify

def classify_feature(name: str, rec: dict) -> tuple[bool, str]:
    """-> (is_gene, basis). Category columns win over name prefixes."""
    if rec.get("is_gene") is not None:
        return rec["is_gene"], "is_gene column"
    cat = rec.get("category")
    if cat:
        c = cat.lower()
        if c in ("gene", "gene expression") or c.endswith("_gene"):
            return True, "category column"
        return False, "category column"
    for prefix in CONTROL_PREFIXES:
        if name.startswith(prefix):
            return False, "name prefix"
    return True, "name prefix"


def control_class(name: str, rec: dict) -> str:
    if rec.get("category"):
        return rec["category"]
    for prefix, label in CONTROL_PREFIXES.items():
        if name.startswith(prefix):
            return label
    return "other_non_gene"


# --------------------------------------------------------------------------- match

def load_library(path=LIBRARY_PATH) -> dict:
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def match_panels(ft: FeatureTable, library: dict) -> dict:
    genes, controls, basis = set(), {}, set()
    for name, rec in ft.features.items():
        is_gene, how = classify_feature(name, rec)
        basis.add(how)
        if is_gene:
            genes.add(name)
        else:
            k = control_class(name, rec)
            controls[k] = controls.get(k, 0) + 1
    if not genes:
        raise ValueError("No gene features were found in this file.")

    ranking = []
    for p in library["panels"]:
        pg = set(p["genes"])
        found = genes & pg
        ranking.append({
            "id": p["id"], "name": p["name"], "species": p["species"], "chemistry": p["chemistry"],
            "design_id": p.get("design_id"), "design_version": p.get("design_version"),
            "source_url": p["source_url"],
            "panel_genes": len(pg),
            "found": len(found),
            "missing": len(pg) - len(found),
            "coverage": len(found) / len(pg),           # fraction of the panel present in the file
            "extra": len(genes - pg),                   # file genes not in this panel
            "jaccard": len(found) / len(genes | pg),
            "_missing_list": sorted(pg - genes),
            "_extra_list": sorted(genes - pg),
        })
    ranking.sort(key=lambda r: (-r["coverage"], -r["jaccard"]))
    best = ranking[0]

    if best["coverage"] >= MATCH_MIN:
        call = "match"
    elif best["coverage"] >= PARTIAL_MIN:
        call = "uncertain"
    else:
        call = "no_match"

    # What the file itself declares (newer XOA writes predesigned_gene / custom_gene)
    declared_custom = sorted(n for n in genes if (ft.features[n].get("category") or "").lower() == "custom_gene")
    has_declared = any((ft.features[n].get("category") or "").lower() in ("custom_gene", "predesigned_gene") for n in genes)

    if call == "no_match":
        custom = sorted(genes)
    else:
        custom = best["_extra_list"]

    return {
        "call": call,
        "best": best,
        "ranking": ranking,
        "n_features": len(ft.features),
        "n_genes": len(genes),
        "genes": sorted(genes),
        "controls": dict(sorted(controls.items())),
        "classification_basis": sorted(basis),
        "custom_genes": custom,
        "missing_genes": best["_missing_list"] if call != "no_match" else [],
        "declared_custom": declared_custom if has_declared else None,
        "source_kind": ft.source_kind,
        "n_rows": ft.n_rows,
        "notes": ft.notes,
        "counts": {n: ft.features[n].get("count") for n in genes},
        "library_version": library.get("library_version"),
    }


def summary_text(res: dict) -> str:
    """Plain-language verdict (Markdown)."""
    b = res["best"]
    pct = f"{b['coverage']:.1%}"
    n_custom = len(res["custom_genes"])
    if res["call"] == "match":
        head = f"### Panel: {b['name']}\n\n"
        body = (f"{b['found']} of the panel's {b['panel_genes']} genes are in the file ({pct}); "
                f"{b['missing']} panel genes are missing.\n\n")
        if n_custom == 0:
            body += "**Custom genes: 0.** Every gene in the file belongs to this pre-designed panel."
        else:
            body += (f"**Custom / add-on genes: {n_custom}.** These genes are in the file but not in the "
                     f"published version of the pre-designed panel ({b['design_id']}), so this looks like the "
                     "pre-designed panel plus a custom add-on. Caveat: 10x has revised some panels, and genes "
                     "that were only in an earlier panel version also show up here, because only the currently "
                     "published version is in the library.")
    elif res["call"] == "uncertain" and n_custom == 0:
        head = f"### Subset of: {b['name']}\n\n"
        body = (f"Every gene in the file belongs to this panel, but only {b['found']} of the panel's "
                f"{b['panel_genes']} genes are present ({pct}); {b['missing']} are missing. "
                "This is not the released panel as published by 10x Genomics. Likely explanations: an earlier "
                "or preview version of the panel, or a file that was filtered after the run. "
                "(Genes with zero detected transcripts are also absent from transcripts.parquet, but that "
                "rarely accounts for more than a few genes.)\n\n**Custom genes: 0.**")
    elif res["call"] == "uncertain":
        head = f"### Uncertain: closest panel is {b['name']}\n\n"
        body = (f"Only {b['found']} of its {b['panel_genes']} genes are in the file ({pct}), below the "
                f"{MATCH_MIN:.0%} needed for a confident call. The file may come from a different panel "
                "version, a filtered export, or an advanced custom design based on this panel. "
                f"{n_custom} genes in the file are not in this panel (listed below); "
                "check the missing-gene list before relying on this.")
    else:
        head = "### No pre-designed panel matches\n\n"
        body = (f"The closest panel ({b['name']}) shares only {b['found']} of its {b['panel_genes']} genes "
                f"with the file ({pct}). This is most likely a fully custom (standalone) panel, or a panel "
                f"that is not in the reference library. All {n_custom} genes in the file are listed below as custom.")
    extra = (f"\n\nFile: {res['source_kind']}, {res['n_features']} distinct features = {res['n_genes']} genes + "
             f"{res['n_features'] - res['n_genes']} control/other features.")
    if res["declared_custom"] is not None:
        dc = set(res["declared_custom"])
        agree = dc == set(res["custom_genes"])
        extra += (f" The file's own `codeword_category` column marks {len(dc)} genes as `custom_gene`"
                  + (", the same set as found here." if agree else
                     ", which differs from the comparison above - see the gene table."))
    return head + body + extra


def custom_table(res: dict) -> list[list]:
    dc = set(res["declared_custom"] or [])
    rows = []
    for g in res["custom_genes"]:
        rows.append([g, res["counts"].get(g), ("yes" if g in dc else "no") if res["declared_custom"] is not None else ""])
    return rows
