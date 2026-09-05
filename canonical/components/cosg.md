# COSG — component reference (self-sufficient)

COSG identifies marker genes (or marker genomic regions) by **cosine specificity**: each gene is
scored for each group by the cosine similarity between its expression vector and the group's
one-hot indicator, then re-scored with the `mu` penalty that down-weights genes also expressed
elsewhere. It is fast (one million cells in under two minutes), works on scRNA-seq, scATAC-seq
and spatial data, and since 1.2.0 carries its own **analytic p-values**. The **same method ships
in two languages**: Python (`cosg`, AnnData **or `.cytome`**) and R (`COSG`, the COSGR repo,
Seurat / SingleCellExperiment). This file assumes nothing about PIASO being installed; COSG is
standalone (**no scanpy dependency** since 1.1.0). Tutorials: https://piaso.org/tutorials/cosg-markers/
(+ `cosg-cytome`, `cosg-batch`, `cosg-gpu`, `cosg-spatial`).

## Install

```bash
pip install cosg                 # Python; numpy/scipy/pandas/numba/anndata only
pip install 'cosg[gpu]'          # + CuPy for device='gpu' | 'auto'
pip install 'cosg[cytome]'       # + cytome for the streaming path (already satisfied if piaso-tools is installed)
```
```r
install.packages("COSG", repos = "https://genecell.r-universe.dev")   # R binaries
# or: conda install -c conda-forge r-cosg      # or: remotes::install_github("genecell/COSGR")
```

Blocks below were executed against **cosg 1.2.0** (Python) and **COSG 1.0.0** (R). Conda's
`bioconda::cosg` lags PyPI — prefer pip.

## Entry points

- **Python:** `import cosg` → `cosg.cosg(data, groupby=..., ...)` — **polymorphic on its first
  argument**: an `AnnData` (in-memory; writes `adata.uns[key_added]`), a **path to a `.cytome`**
  or an **open `CytomeDataset`** (streaming, bounded memory; **returns a dict**, does not close a
  Dataset you passed). `cosg.run_cosg_cytome(...)` is the explicit streaming call. Helpers:
  `indexByGene`, `iqrLogNormalize`, `plotMarkerDotplot`, `plotMarkerDendrogram`, `plotMarkerStream`.
- **R:** `library(COSG)` → `cosg(object, groups = "all", assay = "RNA", slot | layer = "data", ...)`.

## What COSG computes (both languages)

1. One-hot group indicator (groups × cells). 2. Cosine similarity between each gene's expression
vector and each indicator → gene × group score (a gene expressed in exactly one group ≈ 1).
3. `mu` re-scoring: `mu=1` permissive (genes high in this group whether or not high elsewhere);
`mu=100` strict (near-exclusive). It changes the ranking, not the biology. 4. `remove_lowly_expressed`
drops genes detected in fewer than `expressed_pct` of the group's cells (Python also floors at
`expressed_min_num_cells_in_target_group=3`). 5. Top `n_genes_user` per group.

Both expect **normalized, non-negative** values (INFOG layer or log1p) — never a scaled matrix.

## Python block (executed)

Self-sufficient, scanpy-free setup with PIASO's reader is shown; any AnnData with a normalized
layer and a cluster column works.

```python
import pandas as pd
import piaso, cosg                                   # piaso only to read/normalize/cluster; COSG itself needs neither
adata = piaso.pp.read_10x_h5("e18_v3_nuclei.h5")     # 10x h5 with raw counts in .X
piaso.pp.filter_cells(adata, min_counts=500, min_features=250)
piaso.tl.infog(adata, n_top_genes=3000)              # -> layers['infog']
piaso.tl.runSVD(adata, layer="infog", n_components=50, key_added="X_svd")
piaso.tl.neighbors(adata, use_rep="X_svd"); piaso.tl.leiden(adata, resolution=1.0, key_added="leiden")

cosg.cosg(adata, groupby="leiden", key_added="cosg", n_genes_user=25, layer="infog",
          mu=1, expressed_pct=0.1, remove_lowly_expressed=True)
sorted(adata.uns["cosg"].keys())        # ['COSG', 'names', 'params', 'scores']
names = pd.DataFrame(adata.uns["cosg"]["names"])     # one column per cluster, top genes down the rows
marker_sets = {c: list(names[c]) for c in names.columns}
```

`names` / `scores` are structured arrays with one field per group (hence `pd.DataFrame(...)`
gives a group-per-column table). `return_by_group=True` (default) also stores `uns['cosg']['COSG']`,
a wide frame with `names::<group>` / `scores::<group>` columns — the input `indexByGene` expects.
`layer=` is the right habit; do not overwrite `.X`.

Scanpy users: `sc.pp.normalize_total` + `sc.pp.log1p` in `.X` and `sc.tl.leiden` labels work
identically — COSG only needs the normalized matrix and a label column.

### Significance — `calculate_pvalues=True` (Python only)

Analytic p-values on the raw cosine: under the null the cosine is a monotone function of the
gene's within-group sum (a sample sum without replacement) → saddlepoint tail (`pvalue_method=
'spa'`), no permutations, no floor, BH FDR (`pvalue_fdr_method`). Adds `pvals`, `pvals_adj`,
`zscores`, `neg_log10_pvals` to `uns[key_added]`. **P-values are identical across `mu`** — read
FDR as gating association and the score as ranking specificity.

```python
adata.layers["counts"] = adata.X.copy()   # raw UMIs — the null is documented for INTEGER layers
cosg.cosg(adata, groupby="leiden", key_added="cosg_p", n_genes_user=25, layer="counts",
          calculate_pvalues=True)
sorted(adata.uns["cosg_p"].keys())
# ['COSG', 'names', 'neg_log10_pvals', 'params', 'pvals', 'pvals_adj', 'scores', 'zscores']
```

**Double dipping:** when `groupby` came from clustering the same matrix, the p-values are
optimistic (the labels were chosen to separate these cells). Labels from annotation, another
modality or a count split are the valid cases. Many top genes share `pvals = 2.2e-308` — that is
the float floor, not a tie.

### Batches — `batch_key=`

`cosg.cosg(adata, groupby=, batch_key="Sample", batch_cell_number_threshold=5)` computes the
cosine **within each batch and averages**, so a marker must hold in every batch rather than the
largest. Groups with fewer than `batch_cell_number_threshold` cells in a batch are dropped from
that batch's average (a cosine from one cell is 1.0 by construction — raise the threshold when
batches are small). On the cortex tutorial 86 % of top-20 markers survive; the 14 % that move
name the dissociation-sensitive types.

### GPU — `device=`

`cosg.cosg(adata, groupby=, device="gpu")` (or `"auto"`) uses CuPy on both the in-memory and
streaming paths; ~2× in the useful range, slower below ~10,000 cells.

### On a `.cytome` (streaming, returns a dict)

```python
import cytome
ds = cytome.from_anndata(adata, modality="RNA", output="e18.cytome")   # or cytome.open(path) / piaso.data.load_dataset(..., return_type="cytome")
res = cosg.cosg(ds, groupby="leiden", modality="RNA", layer="infog", n_genes_user=25)   # NO key_added — nothing to write to
sorted(res.keys())                                   # ['groups_order', 'names', 'scores']  (+ 'pvals', 'pvals_adj' with calculate_pvalues=True)
markers = pd.DataFrame(res["names"], columns=list(res["groups_order"]))
```

`layer="auto"` (default) normalizes the stored raw counts on the fly; a named layer is read
as-is — pin `layer=` on both sides when comparing with an in-memory run. `output_format` chooses
the shape: `"ndarray"` (default; `names`/`scores` arrays + `groups_order`), `"dict"`
(`scores_dict[(group, gene)]`), `"dense"` (a full feature × group DataFrame — the input for
`iqrLogNormalize`), `"long"`. See `components/cytome.md`.

### Comparing scores across cell types — `iqrLogNormalize`

A raw COSG score is meaningful **within** a group; across groups the scale varies with group size
and exclusivity (112-fold spread of top scores on the cortex tutorial). Run COSG with
`n_genes_user=adata.n_vars`, reshape with `cosg.indexByGene(pd.DataFrame(adata.uns[...]["COSG"]))`
(genes × groups), then `cosg.iqrLogNormalize(scores, q_upper=0.95, q_lower=0.75)` divides each
column by its upper-tail IQR and applies `log1p`. Use it for any shared heatmap, dendrogram or
cross-type threshold. (Whether SCALAR's input should be normalized this way is an open question
on the maintainers' side; `piaso.tl.specificity_matrix` returns the raw matrix.)

### Plots

`cosg.plotMarkerDotplot(adata, groupby=, top_n_genes=2, use_rep="X_umap", key_cosg="cosg",
swap_axes=False, save=<path>)` draws its own layout (returns the `Axes`; `save=` is a path, unlike
scanpy); `cosg.plotMarkerDendrogram(adata, group_by=, use_rep=, top_n_genes=3, cosg_key=)`
clusters groups on their COSG scores; `plotMarkerStream` is a streamgraph. None need scanpy
(`backend="scanpy"` restores the old dot plot).

## R block (executed)

R COSG uses the active `Idents()` (no `groupby`), returns a plain `list(names, scores)`, and does
**not** mutate the object.

```r
library(SeuratObject); library(COSG)
# obj: a Seurat object with a normalized 'data' layer and cell-type Idents()
res <- cosg(obj, groups = "all", assay = "RNA", layer = "data",     # Seurat v5: layer=; v3/v4: slot = "data"
            mu = 1, n_genes_user = 100)
head(res$names)      # one column per identity, top markers down the rows
head(res$scores)
```

`groups = c("0", "2")` restricts to selected identities; `remove_lowly_expressed = TRUE` +
`expressed_pct` as in Python. Seurat v3/v4/v5 and SingleCellExperiment are supported.

## Python ↔ R divergence table

| Aspect | Python `cosg.cosg` (1.2.0) | R `cosg` (COSGR 1.0.0) | Consequence |
|---|---|---|---|
| Data object | `AnnData` (cells × genes) **or `.cytome`** | `Seurat` / SCE (genes × cells) | different orientation / mental model |
| Group labels | `groupby=` column (default `'CellTypes'`) | active `Idents()` — **no `groupby`** | R user must `Idents(obj) <-` first |
| Expression source | `layer=` → `raw.X` (`use_raw`) → `.X` | `assay` + `slot` / `layer` (v5) | R-only params |
| `mu` | `1` | `1` | same |
| `remove_lowly_expressed` | **`True`** (since 1.1.0) | `TRUE` | **same** (they used to differ) |
| `expressed_pct` | `0.1` (+ floor of 3 cells) | `0.1` (no floor) | minor |
| `n_genes_user` | **`50`** | **`100`** | **default differs → different marker counts** |
| Output | writes `adata.uns[key_added]` (`names`, `scores`, `params`, `COSG`; + p-value arrays), returns `None`; on a cytome **returns a dict** | returns `list(names, scores)`; object untouched | different plumbing |
| Python-only | `calculate_pvalues`, `batch_key`, `device='gpu'`, cytome streaming, `calculate_logfoldchanges`, `reference`, `copy`, `output_format`; `indexByGene` / `iqrLogNormalize`; the three plot functions | absent — R exports `cosg()` only | R users lack p-values, batch mode, streaming, plots |

## Python-vs-R disambiguation rule

Infer the language from the objects in the session: `.h5ad` / `AnnData` / `.cytome` / scanpy →
**Python `cosg`**; Seurat object / `.rds` / `library(Seurat)` → **R `COSG`**. Ask when genuinely
ambiguous. An R user who needs p-values or batch-aware markers can write a `.cytome` with the R
`cytome` package (`write_cytome(obj, path)`) and run Python `cosg.cosg(path, groupby=, ...)` on it.

## Citation

Same paper for both implementations:

> Dai M, Pei X, Wang X-J. Accurate and fast cell marker gene identification with COSG.
> *Briefings in Bioinformatics* 23(2):bbab579 (2022). DOI: 10.1093/bib/bbab579

(Both repo READMEs omit the "23(2)" volume/issue; it is restored here. The article ID is bbab579,
not bbac157.)
