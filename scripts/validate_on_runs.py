#!/usr/bin/env python3
"""Check the panel call against what each Xenium run records about itself.

    python scripts/validate_on_runs.py <xenium_outs_dir> [<xenium_outs_dir> ...]

For every output folder: reads transcripts.parquet (read-only), calls the panel,
and prints it next to the panel name / target counts in experiment.xenium.
"""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import panel_core as pc  # noqa: E402

lib = pc.load_library()
for d in sys.argv[1:]:
    d = Path(d)
    truth = {}
    try:
        x = json.loads((d / "experiment.xenium").read_text(encoding="utf-8"))
        truth = {k: x.get(k) for k in ("panel_name", "panel_design_id", "panel_num_targets_predesigned",
                                       "panel_num_targets_custom", "analysis_sw_version")}
    except Exception as e:  # noqa: BLE001
        truth = {"error": str(e)}
    t = time.time()
    try:
        f = d / "transcripts.parquet"
        res = pc.match_panels(pc.read_parquet_features(f), lib)
        b = res["best"]
        out = {"run": d.name, "file_GB": round(f.stat().st_size / 1e9, 2), "seconds": round(time.time() - t),
               "call": res["call"], "best": b["id"], "found": b["found"], "panel_genes": b["panel_genes"],
               "coverage": round(b["coverage"], 4), "n_genes_in_file": res["n_genes"],
               "n_custom": len(res["custom_genes"]), "custom_first10": res["custom_genes"][:10],
               "missing": res["missing_genes"][:10], "basis": res["classification_basis"],
               "declared_custom": None if res["declared_custom"] is None else len(res["declared_custom"]),
               "runner_up": [res["ranking"][1]["id"], round(res["ranking"][1]["coverage"], 3)],
               "experiment.xenium": truth}
    except Exception as e:  # noqa: BLE001
        out = {"run": d.name, "error": f"{type(e).__name__}: {e}", "experiment.xenium": truth}
    print(json.dumps(out), flush=True)
