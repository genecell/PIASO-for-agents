# Ecosystem-wide gotchas

A tight reference of traps that span the PIASO ecosystem. Component-specific detail lives in
each `components/*.md`. Every item here either bit a published tutorial or a hub test.

## Install / import

- **Do not pin `matplotlib<3.9`.** piaso-tools ≥ 1.2.0 and cosg ≥ 1.1.2 import on current
  matplotlib (verified on 3.11.1). Only piaso-tools ≤ 1.1.0 had the removed `cm.get_cmap`
  call; if you meet an old environment that fails on `import piaso`, upgrade (`pip install -U
  piaso-tools`) rather than downgrading matplotlib.
- **`pip install piaso-tools` is the whole install.** It pulls `cosg` and `cytome` as hard
  dependencies; scanpy is **not** installed and **not** needed (`piaso-tools[scanpy]` is an
  interop extra). `pip install scrublet` and `scanpy[leiden]` are unnecessary: `piaso.pp.scrublet`
  and `piaso.tl.leiden` (igraph) are built in.
- **Conda channels lag PyPI.** `bioconda::piaso` was 1.0.3 and `bioconda::cosg` 1.1.3 when
  piaso-tools 1.2.3 / cosg 1.2.0 were on PyPI; laris, cytome, cytorete and emergene have no conda
  package. Prefer `pip`. R packages come from `https://genecell.r-universe.dev` (binaries) or
  conda-forge `r-cosg`.
- **Minimum versions the hub's blocks assume:** piaso-tools 1.2.3 (`predictCellTypeByGDR` silently
  discarded its result on AnnData in 1.2.0–1.2.2), cosg 1.1.2 (matplotlib fix, plain AnnData no
  longer needs cytome), laris 0.13.0 (p-value model), cytome 0.3.0 (`counts` invariant).
- **Emergene pins `annoy < 1.17.0`; PIASO `stitchSpace` (BBKNN) segfaults with newer annoy.**
  `pip install emergene` caps annoy for this reason; pin `annoy==1.16.3` when running
  `stitchSpace` (which also rejects `use_rep` embeddings containing NaN/Inf).

## Deprecated names an agent may meet in old notebooks

| Old (still runs, warns) | Current | Note |
|---|---|---|
| `piaso.tl.runSVDLazy(...)` | `piaso.tl.infog(adata, n_top_genes=)` then `piaso.tl.runSVD(adata, layer="infog")` (or `piaso.tl.infog_svd`) | `runSVDLazy` is a deprecated alias of `infog_svd` |
| `laris.tl.runLARIS(lr_adata, adata=...)` | `laris.tl.runLARIS(lr_data, data, ...)` | `data=` also accepts a cytome Dataset or path |
| `runLARIS(..., n_permutations=, n_repeats=)` | `bg = prepareLRBackground(...)`; `runLARIS(..., background=bg)` | `n_permutations` is legacy-only; passing it with `background=` raises a `FutureWarning`; `n_repeats` is ignored |
| `sc.tl.leiden(adata, key_added="Leiden")` in PIASO tutorials | `piaso.tl.leiden(adata)` → `obs["leiden"]` (lowercase) | old 1.1.0 tutorials used `Leiden`; downstream defaults now assume `leiden` |
| `cosg.run_cosg_cytome_cpu(...)` / `cytome_layer=` | `cosg.cosg(path_or_ds, ...)` or `cosg.run_cosg_cytome(...)` / `layer=` | old names raise `TypeError` with a migration hint |
| `piaso.tl.inferRegulon` / `regulonActivity` / `inferGRN` | `import cytorete; cytorete.inferRegulon(...)` | the `piaso.tl` names are forwarders that raise `ImportError` if cytorete is missing |
| `adata.X = adata.layers["infog"]; cosg.cosg(adata)` | `cosg.cosg(adata, layer="infog")` | every function takes `layer=`; leave `.X` as the author stored it |

## Data / layer contracts

- **INFOG (`piaso.tl.infog`) needs RAW INTEGER counts** (`adata.X` by default, or `layer=`). It
  refuses non-integer input and names the layer to pass (`allow_non_integer=True` overrides).
  Registry datasets may keep counts in a layer — check
  `piaso.data.dataset_info(name)["counts_layer"]` (e.g. `sea_ad_mtg_20k` → `UMIs`,
  `adult_cortex_multiome_rna` → `raw`). Passing the wrong layer to `infog` raises; passing it to
  `score` does **not** — it just resolves worse.
- **`runSVD` defaults to `adata.X`.** Always pass `layer="infog"`, otherwise the SVD runs on raw
  counts and silently ignores the normalization.
- **`piaso.tl.score` defaults to `layer='infog'`** (normalized) and errors if that layer is absent
  — run `infog` first. Single gene set → writes `adata.obs[key_added]` + `adata.uns[key_added]`
  (a DataFrame with 11 score / p-value columns when `compute_pvalues=True`) and returns `None`;
  a dict / DataFrame / list-of-lists → **returns a tuple** `(score_matrix, names, pvals)`.
- **COSG reads normalized values** (`layer="infog"` or log1p in `.X`; non-negative). **COSG
  p-values (`calculate_pvalues=True`) need an integer counts layer** — the saddlepoint null is
  built on raw sums, so pass `layer="counts"` (or wherever the UMIs are) for the p-value run.
  P-values are identical across `mu`; they are **optimistic when `groupby` came from clustering
  the same matrix** (double dipping) — labels from annotation, another modality or a count split
  are the valid cases.
- **`runGDR` defaults changed**: `n_gene=20`, `mu=10.0`, PIASO's own scorer. With `batch_key`,
  pass `infog_layer=<counts layer>` when `.X` is not raw counts (GDR re-runs INFOG per batch;
  pointing it at a scaled `.X` produces NaN and a failed SVD). Output width = number of marker
  groups, not `n_svd_dims`.
- **`getMarkers(..., as_dict=True)` returns a TUPLE** `(markers_df, marker_sets)`. Assigning it
  to one name hands the tuple to `predictCellTypeByMarker` — a mistake that reached three
  published tutorials. Cell-type names must match the study's vocabulary exactly
  (`"319 Astro-TE NN"`, not `"Astro-TE"`); use `list_cell_types=True` first.
- **`predictCellTypeByMarker(use_rep=)` defaults to `X_gdr`.** If your pipeline produced only
  `X_svd`, pass `use_rep="X_svd"`; smoothing needs that embedding to exist. Its labels are the
  **keys** of `marker_gene_set` — cluster IDs if that is what you keyed by.
- **Species-cased QC prefixes.** `calculateCellMetrics(prefix_vars={"mt": "mt-"})` for mouse,
  `"MT-"` for human (`Rps`/`Rpl` vs `RPS`/`RPL`). Matching nothing does not raise — it yields a
  column of zeros that reads like a clean sample.
- **Re-run INFOG → SVD → neighbors → Leiden after removing cells/clusters.** The gene selection
  and SVD were fitted with the removed cells present.

## cytome

- **`{modality}_counts` holds raw integers or does not exist** (cytome ≥ 0.3.0). A non-integer
  `adata.X` is stored as `{modality}_data`; `cosg.cosg(path, layer="auto")` reads the recorded
  source layer. Older files written from normalized `.X` had `log1p` applied twice by
  `layer='auto'` before cosg 1.1.3 — re-run those.
- **`layer="auto"` normalizes; a named layer is read as-is.** Comparing a cytome COSG run with an
  in-memory one and seeing different rankings is almost always this — pin `layer=` on both.
- **Embeddings are renamed on conversion**: `obsm['X_umap']` → `RNA_umap`, `obsm['spatial']` →
  `RNA_spatial` (the `X_` / `obsm` tokens are dropped; `to_anndata` restores them). Files written
  before 0.2.6 use `RNA_obsm_X_umap`; both are readable.
- **One writer per file.** Use the Dataset that `from_10x_h5` / `from_anndata` returns; calling
  `cytome.open(path)` on the next line leaves two writers and the next write fails with
  `database is locked`. Call `ds.flush()` after writing columns/embeddings and `ds.close()` when done.
- **`cosg.cosg(ds)` takes no `key_added`** (there is no `.uns`); it returns a dict whose columns
  are in `groups_order` — re-index any per-cluster table by that order before zipping.
- **Subsetting a cytome is a copy**: `piaso.pp.filter_cells(ds, mask=, inplace=False, output=...)`
  writes a new file; `cell_mask=` on plots/functions selects without writing.
- **`cytome.from_10x_h5` needs `h5py`**, which `pip install cytome` alone does not declare (0.3.1) — install
  `h5py` (or `cytome[anndata]`) in a cytome-only environment.
- `save_layer=True` on `infog` persists the normalized matrix on the file (costs disk, saves
  recomputation); `save_layer=False` normalizes on the fly.

## LARIS (spatial ligand–receptor)

- **Coordinates key is `X_spatial`** (`use_rep_spatial="X_spatial"`), not scanpy's `spatial` —
  copy or pass the key accordingly.
- **P-values changed in 0.12 and 0.13.** Exact p-values need `bg = prepareLRBackground(adata,
  lr_df, use_rep_spatial=)` (the expensive, reusable step) and `runLARIS(..., background=bg)`;
  the floor is 1/10,001 at defaults, deterministic, no seed. Results from < 0.13 involving highly
  expressed genes are not comparable.
- **Defaults changed in 0.10**: `mu=0.25`, `sigma='adaptive'`, neighbourhood sizes 20,
  `spatial_weight=3.0`. Results differ from < 0.10 by design.
- Absent LR genes are now **dropped with a warning** (`unmatched='drop'`); pass
  `unmatched='error'` for a custom database. The `groupby` column must have no missing values.
- `runLARIS` returns a **tuple** `(laris_lr, res)` with `by_celltype=True` (default).
- Group by clusters, not by hundreds of fine types: LARIS tests every sender–receiver pair.

## cytorete

- Needs reference data the package does not ship: `piaso.data.fetch_2bit(genome)` (~800 MB),
  `piaso.data.fetch_genome(genome)` (TSS BED) and `piaso.data.fetch_jaspar()`; `hg38` / `mm10`.
- Only the **RNA regulon chain** ships (`inferRegulon` → `regulonActivity` / `regulonSpecificity`);
  `inferGRN`, `inferTFActivity`, `pp.build_peak_cistrome` raise `ImportError` at call time.
- Per-cell regulon p-values floor at `1/(n_ctrl_set+1)` (−log10 p saturates at 3.0 with the
  default 1,000 control sets) — saturated panels mean "as significant as the test can report".

## PIASOmarkerDB

- **PIASOmarkerDB is a REMOTE REST API — it needs internet.** `piaso.tl.getMarkers` /
  `queryPIASOmarkerDB` / `analyzeMarkers` issue HTTP calls to
  `https://piaso.org/piasomarkerdb/api/v1/`; caches under `~/.piaso/markers`. Python-only client;
  R / other agents use the `piaso-mcp` server or the REST API directly.
