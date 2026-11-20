"""Xenium Panel Identifier - Gradio app (SciLifeLab Serve: Gradio app type, port 7860)."""
from __future__ import annotations

import csv
import os
import re
import tempfile
from pathlib import Path

import gradio as gr

import panel_core as pc

LIBRARY = pc.load_library()
MAX_UPLOAD_MB = os.environ.get("MAX_UPLOAD_MB")  # e.g. "100" on SciLifeLab Serve; unset = no app-side limit

INTRO = f"""
# Xenium Panel Identifier

Upload a 10x Genomics Xenium `transcripts.parquet` (for example a `GSM…_transcripts.parquet` from GEO) to find out
which pre-designed Xenium gene panel the run used and which genes were added as custom targets.

Reference library: {len(LIBRARY['panels'])} pre-designed panels parsed from the official 10x Genomics panel files
(library version {LIBRARY['library_version']}; see the last tab).

**File too large?** On SciLifeLab Serve uploads are limited to 100 MB. For bigger files, either upload the sample's
`features.tsv.gz` / `gene_panel.json` instead (GEO often provides them), or run `extract_features.py` from the
source repository on your computer and upload the small `.features.tsv` it writes.
"""

ACCEPTED = [".parquet", ".pq", ".tsv", ".csv", ".txt", ".gz", ".json"]


def _write_tsv(path, header, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)
    return str(path)


def analyze(file_path, progress=gr.Progress()):
    if not file_path:
        raise gr.Error("Please upload a file first.")
    try:
        progress(0, desc="Reading feature names")
        ft = pc.read_any(file_path, progress=lambda f: progress(f, desc="Reading feature names"))
        res = pc.match_panels(ft, LIBRARY)
    except gr.Error:
        raise
    except Exception as e:  # noqa: BLE001 - show the reason to the user
        raise gr.Error(f"Could not analyse this file: {e}") from e

    ranking = [[r["name"], r["species"], r["chemistry"], r["panel_genes"], r["found"], r["missing"],
                f"{r['coverage']:.1%}", r["extra"]] for r in res["ranking"]]
    custom = pc.custom_table(res)

    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", Path(file_path).name)
    stem = re.sub(r"\.(parquet|pq|tsv|csv|txt|json|gz)$", "", stem, flags=re.I) or "xenium"
    out = Path(tempfile.mkdtemp(prefix="xpi_"))
    f_custom = _write_tsv(out / f"{stem}.custom_genes.tsv",
                          ["gene", "n_transcripts", "marked_custom_gene_in_file"], custom)
    f_missing = _write_tsv(out / f"{stem}.missing_panel_genes.tsv", ["gene"], [[g] for g in res["missing_genes"]])
    f_rank = _write_tsv(out / f"{stem}.panel_ranking.tsv",
                        ["panel", "species", "chemistry", "panel_genes", "found_in_file", "missing_from_file",
                         "panel_coverage", "file_genes_not_in_panel"], ranking)

    ctrl = ", ".join(f"{k}: {v}" for k, v in res["controls"].items()) or "none"
    details = (f"**Controls and other non-gene features** ({res['n_features'] - res['n_genes']}): {ctrl}  \n"
               f"**Genes separated from controls using:** {', '.join(res['classification_basis'])}  \n"
               + (f"**Transcript rows in file:** {res['n_rows']:,}  \n" if res["n_rows"] else "")
               + "".join(f"**Note:** {n}  \n" for n in res["notes"]))

    custom_head = f"### Custom / add-on genes: {len(custom)}"
    missing_txt = (f"### Panel genes missing from the file: {len(res['missing_genes'])}\n\n"
                   + (", ".join(res["missing_genes"]) if res["missing_genes"] else "None."))
    return (pc.summary_text(res), ranking, custom_head, custom, ", ".join(res["custom_genes"]),
            [f_custom, f_missing, f_rank], missing_txt, details)


def library_rows():
    return [[p["name"], p["species"], p["chemistry"], p["n_genes"], f"{p.get('design_id')} v{p.get('design_version')}",
             p["retrieved"], p["source_url"]] for p in LIBRARY["panels"]]


with gr.Blocks(title="Xenium Panel Identifier") as demo:
    gr.Markdown(INTRO)
    with gr.Row():
        inp = gr.File(label="transcripts.parquet, or a feature list (.features.tsv, features.tsv.gz, gene_panel.json)",
                      file_types=ACCEPTED, type="filepath")
    btn = gr.Button("Identify panel", variant="primary")

    with gr.Tabs():
        with gr.Tab("Result"):
            verdict = gr.Markdown("Upload a file and press **Identify panel**.")
            gr.Markdown("#### All panels, best match first\n*Panel coverage* = share of the panel's genes found in "
                        "the file. *File genes not in panel* = what would count as custom if that panel were the base.")
            rank = gr.Dataframe(headers=["Panel", "Species", "Chemistry", "Panel genes", "Found in file",
                                         "Missing from file", "Panel coverage", "File genes not in panel"],
                                interactive=False, wrap=True)
            details = gr.Markdown()
        with gr.Tab("Custom genes"):
            custom_head = gr.Markdown("### Custom / add-on genes")
            custom_df = gr.Dataframe(headers=["Gene", "Transcripts in file", "Marked custom_gene in file"],
                                     interactive=False)
            custom_txt = gr.Textbox(label="Comma-separated list", lines=4, interactive=False)
            files = gr.File(label="Downloads (TSV): custom genes, missing panel genes, panel ranking",
                            file_count="multiple", interactive=False)
        with gr.Tab("Missing panel genes"):
            missing = gr.Markdown()
        with gr.Tab("Reference library"):
            gr.Markdown(
                f"Library version **{LIBRARY['library_version']}**, built by `scripts/build_panel_library.py` from the "
                "gene_panel JSON files linked on the 10x Genomics pages "
                + " and ".join(f"[{u.rsplit('/', 1)[-1]}]({u})" for u in LIBRARY["index_pages"])
                + ". Panels that could not be retrieved: "
                + (", ".join(p["name"] for p in LIBRARY["not_retrieved"]) or "none")
                + ".\n\nNot covered: fully custom (standalone / advanced) panels, Xenium Protein subpanels, and "
                  "pre-release panels that 10x never listed on those pages (for example the 4,624-gene "
                  "*Xenium Human 5K Preview Panel* used in some early 10x demo datasets). "
                  "A call of *match* needs at least "
                  f"{pc.MATCH_MIN:.0%} of a panel's genes in the file; below {pc.PARTIAL_MIN:.0%} is reported as no match.")
            gr.Dataframe(value=library_rows(), headers=["Panel", "Species", "Chemistry", "Genes", "Design",
                                                        "Retrieved", "Source file"], interactive=False, wrap=True)

    btn.click(analyze, inputs=inp,
              outputs=[verdict, rank, custom_head, custom_df, custom_txt, files, missing, details])

if __name__ == "__main__":
    demo.queue(max_size=20).launch(
        server_name=os.environ.get("GRADIO_SERVER_NAME", "127.0.0.1"),
        server_port=int(os.environ.get("GRADIO_SERVER_PORT", "7860")),
        max_file_size=f"{MAX_UPLOAD_MB}mb" if MAX_UPLOAD_MB else None,
    )
