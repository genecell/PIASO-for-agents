# LARIS — component reference (self-sufficient)

LARIS (**L**igand **A**nd **R**eceptor **I**nteraction in **S**patial
transcriptomics) infers spatially-specific ligand–receptor interactions and
sender→receiver cell-type communication for **spatial** data (Visium, MERFISH,
Xenium, Slide-seq, Stereo-seq, etc.). It keys off a spatial kNN graph and is
meaningless without spatial coordinates. This file assumes nothing about PIASO
being installed; LARIS is a standalone package (it does pull in `cosg` as a
dependency, used internally for cell-type specificity).

## Install

```bash
pip install laris        # also pulls cosg (used internally for cell-type specificity)
```

## Import / entry points

`import laris` (or `import laris as la`). Short submodule aliases:
`laris.tl` = `laris.tools`, `laris.pp`, `laris.pl`, and `laris.datasets`.

Core: `laris.tl.prepareLRInteraction`, `laris.tl.runLARIS`.
Bundled LR database: `laris.datasets.lrDatabase`.

## What it computes

- **`prepareLRInteraction(adata, lr_df, number_nearest_neighbors=10,
  use_rep_spatial='X_spatial')`** — builds a spatial kNN graph over the coordinates,
  diffuses expression across it, and for each LR pair takes the element-wise product
  of the diffused ligand × diffused receptor. Returns a **new AnnData** (cells ×
  LR-pairs) with `.var_names = "ligand::receptor"` and spatial coords carried over.
  Requires spatial coordinates in **`adata.obsm['X_spatial']`** (LARIS's default key,
  NOT scanpy's `'spatial'`); `lr_df` must have `ligand`/`receptor` columns whose gene
  names exist in `adata.var_names`.
- **`runLARIS(lr_adata, adata=..., groupby='CellTypes', ...)`** — Step 1 scores each
  LR pair's **spatial specificity** (cosine of the LR score vs its spatially-diffused
  version, minus a shuffled-graph null). Step 2 (`by_celltype=True`, default) uses
  COSG cell-type specificity of ligand/receptor genes plus spatial neighborhood
  cell-type co-localization to produce per sender→receiver interaction scores with
  permutation p-values (BH-FDR). With `by_celltype=True` (default), `adata` is
  **required** and must carry `.obs[groupby]` cell-type labels.

## Bundled LR database (CellChatDB)

`laris.datasets.lrDatabase(species=...)` loads a bundled, curated **CellChatDB**:
**human = 2951** LR pairs, **mouse = 3105** LR pairs. The returned DataFrame has
`ligand`/`receptor` columns (plus pathway/annotation metadata) and feeds directly
into `prepareLRInteraction`. Pick the species that matches your `adata.var_names`.

## Verified block (verbatim, with both gotchas inline)

Requires spatial coordinates in `adata.obsm['X_spatial']`.

```python
import laris
lrdb = laris.datasets.lrDatabase(species="mouse")            # bundled CellChatDB (3105 mouse pairs)
# GOTCHA 1: filter LR-DB to pairs whose ligand AND receptor are both present in the data
present = set(adata.var_names)
lrdb_f = lrdb[lrdb["ligand"].isin(present) & lrdb["receptor"].isin(present)].copy()
# GOTCHA 2: the groupby column must be a categorical dtype
adata.obs["CellTypes"] = pd.Categorical(adata.obs["Leiden"].astype(str))
lr_adata = laris.tl.prepareLRInteraction(adata, lr_df=lrdb_f, number_nearest_neighbors=10,
                                         use_rep_spatial="X_spatial")
res = laris.tl.runLARIS(lr_adata, adata=adata, groupby="CellTypes",
                        n_permutations=100, n_top_lr=500, calculate_pvalues=True)
# returns DataFrame (or tuple); cols: ligand, receptor, score, Rank
```

## Decision rule — LARIS vs SCALAR

Both answer "which ligand–receptor interactions / cell–cell communications are
happening?" Pick by whether the data is spatial:

- **Spatial coordinates present** (Visium/MERFISH/Xenium; `.obsm['X_spatial']` /
  `.obsm['spatial']`) → **LARIS** (`laris.tl.runLARIS`).
- **Dissociated single-cell, no coordinates** → **SCALAR** (`piaso.tl.runSCALAR`;
  see `components/piaso.md`). SCALAR needs a user-supplied LR-pair list + specificity
  matrix (no bundled DB), whereas LARIS bundles CellChatDB.

## Citation

LARIS is a **preprint** (bioRxiv, 2025) — a separate publication from the *Nature* (2026) paper.

> M. Dai, T. Török, D. Sun, et al. LARIS enables accurate and efficient ligand and
> receptor interaction analysis in spatial transcriptomics. *bioRxiv* (2025).
> DOI: 10.1101/2025.11.26.690796
