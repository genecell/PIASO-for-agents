# PIASO ecosystem — overview (tier-1 router)

This file is the routing layer for the PIASO single-cell omics ecosystem. Read it first to
decide **which component answers a request**, then open the matching `components/*.md` (self-
sufficient per-tool reference) or `workflows/*.md` (a multi-step analysis task). If an agent
can read only one file, it should be this one.

## What the PIASO ecosystem is

PIASO is a family of **five installable packages plus one data repository** for single-cell and
spatial transcriptomics, from the Gord Fishell Lab (HMS / Broad). The packages are: **PIASO**
(`piaso-tools`, Python) — the umbrella toolkit for INFOG normalization, dimensionality reduction
(SVD, marker-gene-guided GDR), gene-set scoring, cell-type prediction, the PIASOmarkerDB client,
and **SCALAR** single-cell ligand–receptor inference; **COSG** (`cosg`, Python) — fast marker-gene
identification by cosine-similarity scoring; **COSGR** (`COSG`, **R**) — the R implementation of
COSG for Seurat / SingleCellExperiment objects; **LARIS** (`laris`, Python) — ligand–receptor
interaction analysis for **spatial** transcriptomics; and **Emergene** (`emergene`, Python) —
individual-cell differential transcriptomics **across conditions**. The data repository
**PIASO-data** hosts tutorial datasets (Zenodo) and genome references. The ecosystem is
**multi-language**: everything is Python except COSGR, which brings marker identification to R
users (Seurat/`.rds` workflows). The common substrate is an AnnData object (`.h5ad`) in Python.

## Task → component routing table

| The user wants to… | Route to | Where |
|---|---|---|
| Find marker genes for clusters (AnnData / scanpy) | **COSG** — `cosg.cosg` | `components/cosg.md` |
| Find marker genes for clusters (Seurat / `.rds` / R) | **COSGR** — `cosg()` in R | `components/cosg.md` |
| Score a gene set / gene-set enrichment per cell | **PIASO** — `piaso.tl.score` | `components/piaso.md` |
| Normalize raw UMI counts (information-content) | **PIASO** — `piaso.tl.infog` | `components/piaso.md` |
| Dimensionality reduction / SVD embedding | **PIASO** — `piaso.tl.runSVDLazy` | `components/piaso.md` |
| Marker-gene-guided DR / batch integration | **PIASO** — `piaso.tl.runGDR` / `stitchSpace` | `components/piaso.md` |
| Annotate cell types from a marker set | **PIASO** — `piaso.tl.predictCellTypeByMarker` | `workflows/marker_based_annotation.md` |
| Reference-based label transfer | **PIASO** — `piaso.tl.predictCellTypeByGDR` | `components/piaso.md` |
| Infer cell types for a gene list from a curated DB | **PIASO + PIASOmarkerDB** — `analyzeMarkers` | `workflows/markerdb_annotation.md` |
| Ligand–receptor / cell–cell communication, **dissociated scRNA-seq** | **PIASO / SCALAR** — `piaso.tl.runSCALAR` | `workflows/ligand_receptor.md` |
| Ligand–receptor / cell–cell communication, **spatial** | **LARIS** — `laris.tl.runLARIS` | `workflows/ligand_receptor.md` |
| Differential expression **across ≥2 conditions** | **Emergene** — `emergene.tl.runEMERGENE` | `components/emergene.md` |
| Marker/variable genes in **one condition** | **Emergene** — `emergene.tl.runMarkG` | `components/emergene.md` |
| Full clustering pipeline (load → markers → GDR) | end-to-end scRNA-seq | `workflows/end_to_end_scrnaseq.md` |

## Cross-component decision rules — the hub's core value

No single package can state these; they exist because two tools can answer the same request. Each
rule is symmetric — check it in **both directions**.

### 1. Ligand–receptor: SCALAR (single-cell) vs LARIS (spatial)

Key off **spatial coordinates**.
- **Spatial coordinates present** — a Visium / Slide-seq / MERFISH / Xenium / Stereo-seq object,
  coordinates in `adata.obsm['spatial']` or `adata.obsm['X_spatial']`, or the user asks for
  spatially-specific / neighborhood interactions → **LARIS** (`laris.tl.runLARIS`). LARIS bundles
  CellChatDB (human 2951 / mouse 3105 pairs) via `laris.datasets.lrDatabase`.
- **No coordinates — dissociated single-cell RNA-seq** (only expression + cell-type labels) →
  **PIASO SCALAR** (`piaso.tl.runSCALAR`). SCALAR takes a **user-supplied** LR-pair list and a
  specificity matrix (no bundled DB).
- Tie-breaker: coordinates in `.obsm` → LARIS; no coordinates → SCALAR. LARIS's method keys off a
  spatial kNN graph and is meaningless without coordinates; routing a dissociated dataset to LARIS
  (or a spatial one to SCALAR) is the failure mode to avoid.

### 2. COSG: Python (`cosg`) vs R (COSGR)

Key off the **object type in session**.
- **AnnData / `.h5ad` / scanpy** context → Python `cosg.cosg(adata, groupby=...)`.
- **Seurat object / `.rds` / `library(Seurat)`** context → R `COSG::cosg(object, ...)` (uses the
  active `Idents()`, no `groupby` argument).
- Same method, but **defaults diverge**: `remove_lowly_expressed` is `False` (Python) vs `TRUE`
  (R); `n_genes_user` is `50` (Python) vs `100` (R). A default call therefore returns different
  gene sets across languages — state this when a user compares results. Ask when the object type
  is genuinely ambiguous; answering in the wrong language is worse than not answering.

### 3. Emergene: `runEMERGENE` (≥2 conditions) vs `runMarkG` (1 condition)

Key off the **number of experimental conditions**.
- **≥2 conditions to contrast** (disease vs control, stages) → `emergene.tl.runEMERGENE` (requires
  `condition_key`; uses BBKNN cross-condition diffusion).
- **A single condition / just marker or spatially-variable genes** → `emergene.tl.runMarkG`.
- `runEMERGENE` itself warns and points to `runMarkG` when it detects only one condition.

## Every component installs independently

Each package installs and runs on its own — a COSG-only, LARIS-only, or Emergene-only user is
**first-class**, and each `components/*.md` is written to assume nothing else is installed.

```bash
pip install piaso-tools "matplotlib<3.9"   # PIASO (auto-installs cosg as a hard dependency)
pip install cosg                            # COSG (Python) alone
pip install laris                           # LARIS (pulls cosg)
pip install emergene                        # Emergene (pins annoy<1.17.0)
```
```r
remotes::install_github("genecell/COSGR")   # COSGR (R) — not on CRAN
```

One dependency reality to know: **`piaso-tools` and `laris` both hard-depend on `cosg`** and call
`cosg.cosg(...)` internally, so installing either auto-installs COSG. COSG is **not** re-exported
under the `piaso.*` namespace, though — you still `import cosg` (there is no `piaso.cosg`).
`emergene` and `laris` do **not** depend on `piaso-tools`.

**Mandatory pin:** PIASO 1.1.0 fails to `import` under **matplotlib ≥ 3.9** (its plotting module
calls the removed `matplotlib.cm.get_cmap` at import time). Every PIASO install line and code block
must pin `"matplotlib<3.9"` until fixed upstream. Convention across this hub: use the **short**
submodule aliases (`piaso.tl` / `piaso.pp` / `piaso.pl`; `laris.tl`; `emergene.tl`), which are
runtime-identical to the long `.tools` / `.preprocessing` / `.plotting` forms.

## How PIASO composes with scanpy and scvi-tools

PIASO is **AnnData-native** and designed to slot into an existing scanpy pipeline: it reads and
writes standard `adata.X` / `.layers` / `.obsm` / `.obs` / `.uns` fields, so `sc.pp.neighbors`,
`sc.tl.leiden`, and `sc.pl.*` interoperate directly (the end-to-end workflow chains scanpy
clustering between PIASO steps). PIASO is **complementary to scvi-tools, never a competitor**: it
adds fast, specificity-based marker identification via COSG, INFOG information-content
normalization, marker-gene-guided dimensionality reduction (GDR), and single-cell ligand–receptor
inference (SCALAR) — the last of which has **no scvi-tools equivalent**. A user running scVI for
probabilistic integration can still use COSG for markers, `piaso.tl.score` for gene-set scoring,
and SCALAR/LARIS for communication analysis. Choose per task, not per framework.

## Citations (cite the component actually used)

- **PIASO / Emergene** — Wu, S.J., Dai, M. et al. *Pyramidal neurons proportionately alter cortical
  interneuron subtypes.* Nature (2026). DOI 10.1038/s41586-025-09996-8.
- **COSG / COSGR** — Dai M, Pei X, Wang X-J. *Accurate and fast cell marker gene identification with
  COSG.* Briefings in Bioinformatics 23(2):bbab579 (2022). DOI 10.1093/bib/bbab579.
- **LARIS** — M. Dai, T. Török, D. Sun, et al. *LARIS enables accurate and efficient ligand and
  receptor interaction analysis in spatial transcriptomics.* bioRxiv (2025). DOI
  10.1101/2025.11.26.690796 *(preprint — a separate publication from the Nature (2026) paper)*.

All five packages are BSD-3-Clause. PIASO-data tutorial datasets are redistributed under CC BY 4.0
with attribution to original sources (see `data.md`).
