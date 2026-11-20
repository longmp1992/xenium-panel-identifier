#!/usr/bin/env python3
"""Make a small feature list from a large Xenium transcripts.parquet.

SciLifeLab Serve limits uploads to 100 MB, and transcripts.parquet files are
often several GB. Run this on your own computer, then upload the small TSV it
writes (a few hundred kB at most) to the Xenium Panel Identifier app.

    pip install pyarrow
    python extract_features.py GSM1234567_transcripts.parquet
    # -> GSM1234567_transcripts.features.tsv

Only the feature_name column (plus codeword_category / is_gene when present) is
read, one row group at a time. The input file is never modified.
"""
import argparse
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq


def extract(path):
    with open(path, "rb") as fh:
        pf = pq.ParquetFile(fh)
        names = pf.schema_arrow.names
        if "feature_name" not in names:
            sys.exit(f"No 'feature_name' column in {path}; columns: {', '.join(names)}")
        cols = ["feature_name"] + [c for c in ("codeword_category", "is_gene") if c in names]
        agg = {}
        n = pf.metadata.num_row_groups
        for i in range(n):
            tbl = pf.read_row_group(i, columns=cols)
            arrays = []
            for c in cols:
                col = tbl[c]
                if pa.types.is_dictionary(col.type):
                    col = pc.cast(col, col.type.value_type)
                if pa.types.is_binary(col.type) or pa.types.is_large_binary(col.type):
                    col = pc.cast(col, pa.string())
                arrays.append(col)
            g = pa.table(arrays, names=cols).group_by(cols).aggregate([([], "count_all")]).to_pydict()
            for j, name in enumerate(g["feature_name"]):
                if name is None:
                    continue
                key = (name, g.get("codeword_category", [None] * (j + 1))[j], g.get("is_gene", [None] * (j + 1))[j])
                agg[key] = agg.get(key, 0) + g["count_all"][j]
            print(f"\rrow group {i + 1}/{n}", end="", file=sys.stderr)
        print(file=sys.stderr)
    return agg


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("parquet", type=Path)
    ap.add_argument("-o", "--out", type=Path, help="output TSV (default: <input>.features.tsv in the current folder)")
    a = ap.parse_args()
    out = a.out or Path.cwd() / (a.parquet.stem + ".features.tsv")
    agg = extract(a.parquet)
    with open(out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("feature_name\tcodeword_category\tis_gene\tn_transcripts\n")
        for (name, cat, isg), cnt in sorted(agg.items()):
            fh.write(f"{name}\t{cat or ''}\t{'' if isg is None else str(bool(isg)).lower()}\t{cnt}\n")
    print(f"wrote {out} ({len(agg)} features)")


if __name__ == "__main__":
    main()
