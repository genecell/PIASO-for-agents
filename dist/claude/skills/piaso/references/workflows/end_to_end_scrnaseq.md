# Workflow — end-to-end scRNA-seq (load → clusters → markers → GDR)

A full single-cell RNA-seq pipeline built from PIASO + COSG + scanpy: load a 10x dataset, QC-filter,
INFOG-normalize, reduce dimensions, cluster, find markers, and build a marker-gene-guided embedding.
Every block below is executed and passing on the `e18_v3_nuclei` fixture (see `data.md`).

## Install

```bash
pip install piaso-tools cosg igraph leidenalg "matplotlib<3.9"
```

`piaso-tools` pulls scanpy and (as a hard dependency) `cosg`; `igraph` + `leidenalg` are needed by
the scanpy Leiden step. The `"matplotlib<3.9"` pin is mandatory — PIASO 1.1.0 fails to import
without it.

## State that flows between steps

Each step reads specific AnnData fields and writes new ones; the next step depends on them. Track:
`adata.layers["counts"]` (raw UMIs), `adata.layers["infog"]` (normalized), `adata.obsm["X_svd"]`
(SVD embedding), `adata.obs["Leiden"]` (cluster labels), `adata.uns["cosg"]` (markers),
`adata.obsm["X_gdr"]` (marker-guided embedding).

## Step 1 — Load + QC

Runs first. Reads the 10x `.h5`; writes filtered `adata` and preserves raw counts in a layer for
INFOG (which needs raw UMIs).

```python
import scanpy as sc
adata = sc.read_10x_h5("e18_v3_nuclei.h5")   # PIASO-data: e18_v3_nuclei (10x h5, ~5k cells)
adata.var_names_make_unique()
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.layers["counts"] = adata.X.copy()      # keep raw UMIs
```
**Out:** filtered `adata`; `adata.layers["counts"]` = raw UMIs.

## Step 2 — INFOG normalization

Runs after QC (needs the raw counts from Step 1). Writes a normalized layer and selects HVGs.

```python
import piaso
piaso.tl.infog(adata, layer=None, n_top_genes=2000, key_added="infog")
# writes adata.layers['infog'] (normalized) + adata.var['highly_variable']
```
**Out:** `adata.layers["infog"]`; `adata.var["highly_variable"]`.

## Step 3 — SVD embedding

Runs after INFOG. Reduces the normalized data to 50 dimensions for the neighbor graph.

```python
piaso.tl.runSVDLazy(adata, layer="infog", infog_layer="counts",
                    n_components=50, n_top_genes=2000, key_added="X_svd", random_state=1927)
# writes adata.obsm['X_svd'] (n_obs x 50)
```
**Out:** `adata.obsm["X_svd"]` (n_obs × 50).

## Step 4 — Neighbors + Leiden clustering

Standard scanpy, on the SVD embedding. Writes cluster labels used by every downstream step.

```python
sc.pp.neighbors(adata, use_rep="X_svd", n_neighbors=15)
sc.tl.leiden(adata, resolution=1.0, key_added="Leiden", flavor="igraph", n_iterations=2, directed=False)
```
**Out:** `adata.obs["Leiden"]` (cluster labels).

## Step 5 — COSG marker genes

Runs after clustering. COSG expects normalized/log values in `.X`, so point `.X` at the INFOG layer
first. Writes marker recarrays per cluster.

```python
import cosg
adata.X = adata.layers["infog"]          # COSG expects normalized/log values in .X (or a layer)
cosg.cosg(adata, groupby="Leiden", key_added="cosg", n_genes_user=30, mu=1.0)
# writes adata.uns['cosg'] = {'names': recarray, 'scores': recarray}, one field per cluster
```
**Out:** `adata.uns["cosg"]` (`names`, `scores` recarrays, one field per cluster).

## Step 6 — GDR (marker-gene-guided DR)

Final step. Builds a cell × cluster embedding by scoring every cell against each cluster's markers
(runs COSG internally). Its width = number of clusters, not 50.

```python
piaso.tl.runGDR(adata, groupby="Leiden", n_gene=30, mu=1.0, layer="infog",
                score_layer="infog", scoring_method="scanpy", key_added="X_gdr")
# writes adata.obsm['X_gdr'] (n_obs x n_clusters)
```
**Out:** `adata.obsm["X_gdr"]` (n_obs × n_clusters); `adata.uns["gdr"]`.

## Where to go next

`adata.obsm["X_gdr"]` can drive `sc.pp.neighbors(use_rep="X_gdr")` for UMAP/re-clustering, and both
`adata.uns["cosg"]` and `adata.obsm["X_gdr"]` feed marker-based annotation —
see `marker_based_annotation.md`.
