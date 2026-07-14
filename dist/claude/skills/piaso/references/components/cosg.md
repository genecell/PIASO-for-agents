# COSG — component reference (self-sufficient)

COSG is a fast, specific marker-gene identification method: it scores each gene
for each cluster by **cosine similarity** between the gene's expression vector and
a one-hot cluster indicator, then applies a second-stage `mu` penalty that
down-weights genes also expressed in other clusters (higher `mu` → higher
specificity). It is a faster, more specific alternative to `sc.tl.rank_genes_groups`.

The **same method ships in two languages**: Python (`cosg`, working on AnnData)
and R (`COSG`, aka the COSGR repo, working on Seurat). The algorithm is identical
(cosine similarity + the same `mu` second-stage formula), but the two
implementations have **materially different defaults and data contracts** — see the
divergence table below. This file assumes nothing about PIASO being installed;
COSG is a standalone package.

## Install

```bash
# Python
pip install cosg
```
```r
# R (GitHub only — not on CRAN; needs proxyC, data.table, SeuratObject)
remotes::install_github("genecell/COSGR")
```

## Import / entry point

- Python: `import cosg` → `cosg.cosg(...)` (no submodule aliases).
- R: `library(COSG)` → `cosg(...)`.

## What COSG computes (both languages)

1. Build a one-hot cluster indicator matrix (cluster × cell).
2. Compute **cosine similarity** between each gene's expression vector and each
   cluster indicator column → a (gene × cluster) score. A gene expressed in
   exactly one cluster's cells scores ≈ 1.
3. Second-stage penalized re-scoring (the COSG `mu` penalty):
   - `mu == 1`: `score = cosine_sim**2 / row_sum(cosine_sim**2) * cosine_sim`
   - else: `score = cosine_sim**2 / ((1-mu)*cosine_sim**2 + mu*row_sum) * cosine_sim`
   Larger `mu` → stronger penalty for genes also expressed elsewhere.
4. Optional `remove_lowly_expressed`: sets score to −1 for genes expressed in too
   few cells of the target group.
5. Select the top `n_genes_user` genes per group.

Both implementations expect **log-normalized** expression values.

## Python block (verified)

COSG expects **log-normalized** values in `.X` (or a named layer) and reads cluster
labels from `adata.obs[groupby]`. A self-sufficient, scanpy-only setup (no PIASO required):

```python
import scanpy as sc
import cosg
adata = sc.read_10x_h5("your_10x.h5")         # any 10x .h5 or .h5ad with raw counts
adata.var_names_make_unique()
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)  # standard log-normalization
sc.pp.log1p(adata)
sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata)
sc.tl.leiden(adata, key_added="Leiden", flavor="igraph", n_iterations=2, directed=False)

# COSG on the log-normalized values now in .X
cosg.cosg(adata, groupby="Leiden", key_added="cosg", n_genes_user=30, mu=1.0)
# writes adata.uns['cosg'] = {'names': recarray, 'scores': recarray}, one field per cluster
```

> Already normalized with PIASO's INFOG? Set `adata.X = adata.layers["infog"]` before calling
> COSG instead of the `normalize_total` / `log1p` step. PIASO is optional — COSG needs only
> `cosg` (which pulls in `scanpy`).

Read the top markers per cluster back out of the structured arrays:

```python
names = adata.uns["cosg"]["names"]
marker_set = {cl: [names[cl][i] for i in range(len(names))] for cl in names.dtype.names}
```

## R block (verified)

R COSG uses the active `Idents()` (there is **no `groupby` argument**), returns a
plain `list(names, scores)`, and does **not** mutate the Seurat object.

```r
library(SeuratObject); library(COSG)
# obj: a Seurat object with a normalized 'data' layer and cell-type Idents()
# COSG (R) uses Idents() (NO groupby arg), returns list(names, scores), does NOT mutate obj
res <- cosg(obj, groups = 'all', assay = 'RNA', slot = 'data', mu = 1, n_genes_user = 100)
res$names$A   # top markers for identity 'A'
```

## Python ↔ R divergence table

Same method, two languages, **different defaults and contracts**. A default call in
each language returns *different gene sets* (because both `n_genes_user` and
`remove_lowly_expressed` defaults flip).

| Aspect | Python `cosg.cosg` | R `cosg` (COSGR) | Consequence |
|---|---|---|---|
| Data object | `AnnData` (cells × genes) | `Seurat` (genes × cells) | different mental model / orientation |
| Cluster labels | explicit `groupby=` arg, reads `adata.obs[groupby]` (default `'CellTypes'`) | **no `groupby`** — uses active `Idents()` | R user must `SetIdent` first |
| Expression source | `layer` → `raw.X` (if `use_raw`) → `adata.X` | `assay`/`slot` (default `assay='RNA'`, `slot='data'`); `layer` (Seurat v5) takes precedence | R-only `assay`/`slot` params |
| `mu` (penalty) | `1` | `1` | **same** |
| `remove_lowly_expressed` | **`False`** | **`TRUE`** | **default flips → different gene sets** |
| `expressed_pct` | `0.1` | `0.1` | same |
| Lowly-expressed floor | `max(n_cells*expressed_pct, expressed_min_num_cells_in_target_group=3)` | `n_cells*expressed_pct` only (no floor) | Python-only absolute floor |
| `n_genes_user` | **`50`** | **`100`** | **default flips → different # markers** |
| Output | writes `adata.uns[key_added or 'cosg']` (structured `names`/`scores`/`params`), returns `None` (or copy if `copy=True`) | returns `list(names=df, scores=df)`; object untouched | different result plumbing |
| Python-only features | `batch_key` per-batch averaging, `calculate_logfoldchanges`, `reference`, `use_raw`, `copy`, `return_by_group`, `key_added`; plotting suite `plotMarkerDotplot`/`plotMarkerDendrogram`/`plotMarkerStream`; helpers `indexByGene`/`iqrLogNormalize` | absent — R exports only `cosg` | R users lack batch mode + plots |
| Version | 1.0.4 | 1.0.0 | R lags |

## Python-vs-R disambiguation rule

Infer the language from the objects in the session:

- An `.h5ad` / `AnnData` / `scanpy` context → **Python `cosg`**.
- A `Seurat` object / `.rds` / `library(Seurat)` context → **R `COSG` (COSGR)**.
- Ask when genuinely ambiguous. Answering in the wrong language is worse than not
  firing.

## Citation

Same paper for both implementations:

> Dai M, Pei X, Wang X-J. Accurate and fast cell marker gene identification with
> COSG. *Briefings in Bioinformatics* 23(2):bbab579 (2022).
> DOI: 10.1093/bib/bbab579

(Both repo READMEs omit the "23(2)" volume/issue; it is restored here.)
