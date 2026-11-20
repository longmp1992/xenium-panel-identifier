# Xenium Panel Identifier

A small Gradio app that tells you which pre-designed 10x Genomics Xenium gene panel a dataset was run with,
and which genes were added as custom targets. Upload a Xenium `transcripts.parquet` (for example a
`GSM…_transcripts.parquet` from GEO); the app reads only the distinct feature names, separates genes from
control features, compares the gene set with a library of the official panels, and reports

- the best-matching panel with overlap statistics, and all other panels ranked, so you can judge the call;
- the custom / add-on genes (count, table, downloadable TSV);
- panel genes that are missing from the file.

If fewer than 50% of the closest panel's genes are present it says "no pre-designed panel matches" (likely a
standalone custom panel) instead of forcing a match. Between 50% and 95% the result is reported as uncertain.

## Files

| File | Purpose |
|---|---|
| `main.py` | Gradio app (port 7860) |
| `panel_core.py` | Reading, gene/control classification, matching |
| `extract_features.py` | Local helper: makes a small feature list from a large `transcripts.parquet` |
| `data/xenium_panel_library.json` | Bundled reference library (generated) |
| `scripts/build_panel_library.py`, `scripts/panel_sources.json` | Rebuild the library from the official 10x files |
| `scripts/validate_on_runs.py` | Compare calls with `experiment.xenium` of local Xenium output folders |
| `Dockerfile`, `start-script.sh`, `requirements.txt` | Container for SciLifeLab Serve |

## Accepted input

- `transcripts.parquet` (any file name). Columns used: `feature_name`, plus `codeword_category` and `is_gene`
  when present (newer Xenium Onboard Analysis). Read one row group at a time; memory use stays small.
- `*.features.tsv` written by `extract_features.py`.
- `features.tsv` / `features.tsv.gz` from `cell_feature_matrix`, or the run's `gene_panel.json`.
- A plain text / CSV list with one gene name per line.

Genes are separated from controls with `is_gene` / `codeword_category` when the file has them, otherwise by
name prefix (`NegControlProbe_`, `NegControlCodeword_`, `UnassignedCodeword_`, `DeprecatedCodeword_`, `BLANK_`,
`Intergenic_Region_`, `antisense_`).

## Reference library

`data/xenium_panel_library.json` (version 2026-11-20) holds 11 panels, each parsed from the official
`gene_panel` JSON linked on the 10x Genomics Xenium Panel Designer pages
([v1 panels](https://www.10xgenomics.com/support/software/xenium-panel-designer/latest/tutorials/pre-designed-xenium-v1),
[Prime 5K panels](https://www.10xgenomics.com/support/software/xenium-panel-designer/latest/tutorials/pre-designed-xenium-prime-5k)).
Parsed gene counts equal the numbers stated on those pages for all 11.

| Panel | Species | Chemistry | Genes | Design |
|---|---|---|---|---|
| Xenium Human Brain Gene Expression Panel | Human | v1 | 266 | hBrain_v1.2 |
| Xenium Human Breast Gene Expression Panel | Human | v1 | 280 | hBreast_v1.1 |
| Xenium Human Colon Gene Expression Panel | Human | v1 | 322 | hColon_v1 |
| Xenium Human Immuno-Oncology Profiling Panel | Human | v1 | 380 | hImmune_v1 |
| Xenium Human Lung Gene Expression Panel | Human | v1 | 289 | hLung_v1.1 |
| Xenium Human Multi-Tissue and Cancer Panel | Human | v1 | 377 | hMulti_v1 |
| Xenium Human Skin Gene Expression Panel | Human | v1 | 260 | hSkin_v1 |
| Xenium Mouse Brain Gene Expression Panel | Mouse | v1 | 247 | mBrain_v1.1 |
| Xenium Mouse Tissue Atlassing Panel | Mouse | v1 | 379 | mMulti_v1 |
| Xenium Prime 5K Human Pan Tissue & Pathways Panel | Human | Prime 5K | 5001 | hAtlas_v1.1 |
| Xenium Prime 5K Mouse Pan Tissue & Pathways Panel | Mouse | Prime 5K | 5006 | mAtlas_v1 |

Not retrieved: none of the panels listed on those two pages. Not in the library at all:

- earlier versions of revised panels (10x publishes only the current file), and pre-release panels such as the
  4,624-gene "Xenium Human 5K Preview Panel";
- Xenium Protein subpanels and any panel 10x does not list on the two pages above.

Rebuild (standard library only):

```
python scripts/build_panel_library.py
```

The script downloads each URL in `scripts/panel_sources.json`, checks the gene count against the number on the
10x page, and records failures under `not_retrieved` without inventing gene lists. When 10x releases a new panel,
copy its JSON link from the index page into `panel_sources.json` and re-run. `--discover` tries to find new links
automatically, but www.10xgenomics.com answers scripts with HTTP 429, so in practice check the pages in a browser.

## Known limits (from testing on real runs)

Tested on 8 real Xenium output folders (7 v1 runs, 1 Prime 5K run); see `test_outputs/validation.jsonl`.

- The base panel was identified correctly in all 7 v1 runs (100% of the panel genes present each time).
- Custom-gene counts equal `experiment.xenium` in 3 add-on runs (100, 97, 100 genes).
- In 4 runs made with an earlier panel version the count is higher than the true add-on, because genes that
  10x later removed from the panel are not in the published file and are therefore listed as custom: skin base
  (22 listed, 0 real), skin add-on (122 vs 100), mouse brain (1 vs 0), lung preview (103 vs 100).
  The app states this caveat next to the count.
- The Prime 5K lymph node demo dataset was run with the 4,624-gene Preview panel; the app reports
  "subset of Xenium Prime 5K Human" (92.5% coverage, 0 custom) rather than a match.
- Genes with zero detected transcripts do not appear in `transcripts.parquet`, so a few panel genes can be
  missing in small or low-quality samples.

## Run locally

```
pip install -r requirements.txt
python main.py
```

Then open http://localhost:7860. No upload limit applies unless `MAX_UPLOAD_MB` is set.

With Docker:

```
docker build --platform linux/amd64 -t xenium-panel-identifier:0.1.0 .
docker run --rm -p 7860:7860 xenium-panel-identifier:0.1.0
```

The image sets `MAX_UPLOAD_MB=100` to mirror Serve; add `-e MAX_UPLOAD_MB=` to lift it locally.

The image is built by GitHub Actions (`.github/workflows/docker.yml`, run manually with a tag), which also checks
that the container starts, serves the page and runs as user id 1000, then pushes to
`ghcr.io/longmp1992/xenium-panel-identifier:<tag>`.

Deployed at https://xenium-panel-identifier.serve.scilifelab.se (SciLifeLab Serve, image tag 0.1.0).

## Files larger than 100 MB

SciLifeLab Serve limits uploads to 100 MB, and `transcripts.parquet` is usually 0.2 to 3 GB. Options:

1. Run the helper on your own computer and upload the small TSV it writes (about 0.5 MB for a 3 GB 5K file,
   which took about 3 minutes from a network drive):
   ```
   pip install pyarrow
   python extract_features.py GSM1234567_transcripts.parquet
   ```
2. Upload the sample's `features.tsv.gz` or `gene_panel.json` instead, if GEO provides it.
3. Run the app locally (above), where there is no limit.
4. Ask serve@scilifelab.se for a higher upload limit for the app.

## Publish on SciLifeLab Serve

Requirements taken from the Serve documentation (v5.4.0-beta, read 2026-11-20,
https://serve.scilifelab.se/docs/application-hosting/gradio/): Docker image, public registry, container user id
1000, Gradio on port 7860, default 2 vCPU / 4 GB RAM, 100 MB upload limit, source code must be public, app must
eventually be public. The Dockerfile here follows those.

Steps that need your own accounts:

1. Put this folder in a public Git repository (Serve requires the code to be public).
2. Build and push the image to a public registry with a unique tag (do not reuse `latest`):
   ```
   docker build --platform linux/amd64 -t <dockerhub-user>/xenium-panel-identifier:0.1.0 .
   docker push <dockerhub-user>/xenium-panel-identifier:0.1.0
   ```
   (GitHub Container Registry works too; make the package public.)
3. Log in at https://serve.scilifelab.se with your Swedish university / SciLifeLab e-mail.
4. Projects, then create a project with the "Default project template".
5. On the "Gradio app" card click Create and fill in:
   - Subdomain: for example `xenium-panel-identifier`
   - Hardware: default (2 vCPU, 4 GB) is enough
   - Port: `7860`
   - Image: `<dockerhub-user>/xenium-panel-identifier:0.1.0`
   - Title, description, keywords, source code URL
   - Visibility: Private or Link while testing, Public when ready
6. Click Publish; the app is usually reachable a few minutes after the status shows "Running".
7. To update: push a new tag, then Settings of the app, change the Image field, Update.

10x Genomics panel files are redistributed here as parsed gene lists with their source URLs; check the 10x terms
of use before making the repository public.
