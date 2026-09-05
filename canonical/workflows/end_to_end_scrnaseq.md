# Workflow — end-to-end scRNA-seq (load → QC → clusters → markers → cell types → GDR)

A full single-cell RNA-seq pipeline in **PIASO + COSG only — no scanpy**: read a 10x dataset,
QC both tails, call doublets, INFOG-normalize, embed, cluster, find markers, check for artefact
clusters, name the cell types, and rebuild the embedding from the cell types with GDR. Every
block below is executed and passing on the `e18_v3_nuclei` fixture (see `data.md`); the same
steps run unchanged on a `.cytome` (see `streaming_large_data.md`). It mirrors the executed
piaso.org tutorials *Mouse brain scRNA-seq end to end* and *Human PBMC scRNA-seq end to end*.

## Install

```bash
pip install piaso-tools      # installs cosg and cytome too; no matplotlib pin, no scanpy
```

## State that flows between steps

`adata.X` (raw UMI counts, unchanged throughout) → `obs[n_counts, n_genes, pct_counts_*,
scrublet_score, is_doublet]` → `layers["infog"]` + `var["highly_variable"]` → `obsm["X_svd"]` →
`obsp` graph → `obs["leiden"]` → `obsm["X_umap"]` → `uns["cosg"]` → `obs["CellTypes"]` →
`obsm["X_gdr"]`. Every function reads `layer=` / `use_rep=` explicitly; `.X` is never overwritten.

## Step 1 — Load

```python
import numpy as np, pandas as pd
import piaso, cosg
piaso.settings.set_figure_params(style="cell")           # one house style for every figure

adata = piaso.pp.read_10x_h5("e18_v3_nuclei.h5")         # Cell Ranger h5; raw UMI counts in .X, ~6k mouse brain nuclei
# adata = piaso.data.load_dataset("e18_v3_nuclei")       # same file via the registry (downloads + md5-verifies)
# adata = piaso.pp.read_10x("path/filtered_feature_bc_matrix/")   # MTX directory
```
**Out:** `adata` (5973 × 32285), `var[gene_ids, gene_symbols, feature_types, genome]`.

## Step 2 — QC, both tails

```python
piaso.pp.calculateCellMetrics(adata, prefix_vars={"mt": "mt-", "ribo": ["Rps", "Rpl"]})   # mouse; human: "MT-", ["RPS","RPL"]
adata.obs[["n_counts", "n_genes", "pct_counts_mt", "pct_counts_ribo"]].describe()
piaso.pl.scatter(adata, x="n_counts", y="n_genes", color="pct_counts_mt", logx=True, logy=True, marginals=True)

hi_c = np.percentile(adata.obs["n_counts"], 99)          # too MANY counts is a signal too (doublets)
keep = ((adata.obs["n_counts"] >= 500) & (adata.obs["n_counts"] <= hi_c)
        & (adata.obs["n_genes"] >= 250) & (adata.obs["pct_counts_mt"] <= 5.0))
adata = adata[keep].copy()
```
**Out:** filtered `adata`. Read the distribution before fixing a mitochondrial cut: nuclei sit near
0 % (the cut does nothing), whole cells near 10 %+ (a fixed 10 % would remove healthy cells).
Prefixes are case-sensitive and match nothing silently.

## Step 3 — Doublets

```python
piaso.pp.scrublet(adata, expected_doublet_rate=0.06, random_state=0)   # PIASO's own Scrublet; library_key= for several libraries
adata.obs["doublet_call"] = adata.obs["is_doublet"].astype(bool) | (adata.obs["scrublet_score"] > 0.3)
adata = adata[~adata.obs["doublet_call"]].copy()
piaso.pp.filter_cells(adata, min_counts=500, min_features=250)
```
**Out:** `obs["scrublet_score", "is_doublet", "doublet_call"]`; doublets removed. The automatic
threshold is conservative (it needs a dip in the score distribution); `0.3` is a convention to
read off the scatter and the per-cluster violins below.

## Step 4 — INFOG normalization + informative genes

```python
piaso.tl.infog(adata, n_top_genes=3000)                   # reads RAW counts from .X (or layer=); writes layers['infog'], var['highly_variable']
```
**Out:** `adata.layers["infog"]`, `adata.var["infog_var"]`, `adata.var["highly_variable"]`.

## Step 5 — SVD, neighbours, Leiden, UMAP

```python
piaso.tl.runSVD(adata, layer="infog", n_components=50, key_added="X_svd")   # pass layer= — default is .X (raw counts)
piaso.tl.neighbors(adata, use_rep="X_svd", n_neighbors=15)
piaso.tl.leiden(adata, resolution=1.0, key_added="leiden")                  # igraph; obs['leiden'] lowercase
piaso.tl.umap(adata, use_rep="X_svd")
piaso.pl.embedding(adata, basis="X_umap", color="leiden", legend_loc="both")
```
**Out:** `obsm["X_svd"]` (n × 50), `obs["leiden"]` (22 clusters on the fixture), `obsm["X_umap"]`.

## Step 6 — COSG markers, and which clusters are artefacts

```python
cosg.cosg(adata, groupby="leiden", key_added="cosg", n_genes_user=25, layer="infog")
names = pd.DataFrame(adata.uns["cosg"]["names"]); scores = pd.DataFrame(adata.uns["cosg"]["scores"])

top3 = []
for c in adata.obs["leiden"].cat.categories:
    top3 += [g for g in names[c][:3] if g not in top3]
piaso.pl.dotplot(adata, top3, groupby="leiden", standard_scale="var")          # a good cluster's markers are darker AND larger on the diagonal
piaso.pl.plot_features_violin(adata, ["n_genes", "n_counts", "pct_counts_mt", "scrublet_score"], groupby="leiden")

qc = adata.obs.groupby("leiden", observed=True).agg(
    n_cells=("n_genes", "size"), median_genes=("n_genes", "median"),
    median_scrublet=("scrublet_score", "median"),
    frac_high_scrublet=("scrublet_score", lambda v: float((v > 0.2).mean())))
qc["top_cosg"] = [float(scores[c].iloc[0]) for c in qc.index]                 # how specific the best marker is at all
qc["mt_in_top10"] = [sum(g.startswith("mt-") for g in names[c].head(10)) for c in qc.index]   # are the markers themselves technical?
drop = qc.index[(qc["median_scrublet"] > 0.15) | (qc["frac_high_scrublet"] > 0.5) | (qc["mt_in_top10"] >= 3)]
adata = adata[~adata.obs["leiden"].isin(drop)].copy()
```
**Out:** `uns["cosg"]` (`names`, `scores`, `params`, `COSG`); artefact clusters removed. A cluster's
own markers are its most informative QC metric: a top COSG score near 0.02 with half its cells
above the doublet threshold has no identity (doublets); five `mt-` genes in its top ten is dying
cells. Do **not** drop clusters on gene count alone — endothelium and glia carry less RNA and
score the *highest* specificity.

## Step 7 — Re-run the embedding on the cells that remain

```python
piaso.tl.infog(adata, n_top_genes=3000)
piaso.tl.runSVD(adata, layer="infog", n_components=50, key_added="X_svd")
piaso.tl.neighbors(adata, use_rep="X_svd", n_neighbors=15)
piaso.tl.leiden(adata, resolution=1.0, key_added="leiden")
piaso.tl.umap(adata, use_rep="X_svd")
cosg.cosg(adata, groupby="leiden", key_added="cosg", n_genes_user=25, layer="infog")
```
The gene selection, SVD and graph were fitted with the removed cells present; skipping this
analyses clean cells in a space defined partly by the cells you removed.

## Step 8 — Cell types, not cluster numbers

Two marker sources; compare them (they fail differently).

```python
# Route A — PIASOmarkerDB (curated, live API). as_dict=True returns a TUPLE: unpack both.
markers_df, marker_db = piaso.tl.getMarkers(study="AllenWholeMouseBrain_isocortex", as_dict=True)
piaso.tl.predictCellTypeByMarker(adata, marker_gene_set=marker_db, score_layer="infog",
                                 use_rep="X_svd", key_added="CellTypes_db")       # use_rep defaults to X_gdr — pass the embedding you have
adata.obs["CellTypes_db"].value_counts().head()

# Route B — your own clusters' COSG markers -> what PIASOmarkerDB calls them
names = pd.DataFrame(adata.uns["cosg"]["names"])
results, top_hits = piaso.tl.analyzeMarkers({c: list(names[c]) for c in names.columns}, species="Mouse")
adata.obs["CellTypes_cosg"] = adata.obs["leiden"].map(top_hits)                   # top_hits[cluster] is a cell-type string

piaso.pl.embedding(adata, basis="X_umap", color=["CellTypes_db", "CellTypes_cosg"], ncol=1)
```
**Out:** `obs["CellTypes_db"]` (+ `_raw`, `_smoothed`, `_score`, confidence columns,
`obsm["CellTypes_db_score"]`), `obs["CellTypes_cosg"]`. Labels keep the source taxonomy's names
(`"319 Astro-TE NN"`) — traceable to the study. Check the annotation against its own markers with
`piaso.pl.dotplot(adata, <top genes per type>, groupby="CellTypes_db", standard_scale="var")`, and
per-type QC with `gm = piaso.pp.calculateGroupMetrics(adata, groupby="CellTypes_db");
piaso.pl.plotGroupMetrics(gm, data=adata, groupby="CellTypes_db")`. Route C — an annotated
reference dataset of the same tissue: `cosg.cosg(ref, groupby="CellTypes", layer=...)` → marker
dict → `predictCellTypeByMarker`; see `marker_based_annotation.md`.

## Step 9 — GDR: rebuild the embedding from the cell types

```python
piaso.tl.runGDR(adata, groupby="CellTypes_db", layer="infog", key_added="X_gdr")   # n_gene=20, mu=10 defaults; width = n cell types
piaso.tl.neighbors(adata, use_rep="X_gdr", n_neighbors=15, key_added="gdr")
piaso.tl.umap(adata, use_rep="X_gdr", key_added="X_umap_gdr", neighbors_key="gdr")
piaso.pl.embedding(adata, basis="X_umap_gdr", color="CellTypes_db")
adata.write_h5ad("e18_annotated.h5ad")
```
**Out:** `obsm["X_gdr"]`, `uns["gdr"]`, `uns["gdr_reference"]` (so new samples can be projected
with `piaso.tl.projectGDR`), `obsm["X_umap_gdr"]`. GDR's axes are cell-identity axes, so the
embedding separates types rather than directions of maximum variance; with several samples pass
`batch_key="sample", groupby=None` and it integrates by identity.

## Where to go next

- Several samples → `piaso.pp.scrublet(library_key=)`, `runGDR(batch_key=, groupby=None)`; a
  measured batch effect → `piaso.tl.runHarmony` (needs `harmonypy`).
- Gene sets → `piaso.tl.score` (`components/piaso.md`); ligand–receptor → `ligand_receptor.md`;
  regulons → `gene_regulatory_networks.md`; conditions → `components/emergene.md`.
- Too big for RAM, or want the results on disk → the same steps on a cytome: `streaming_large_data.md`.
- Using scanpy alongside: the object is a plain AnnData — `sc.pp.neighbors(adata, use_rep="X_svd")`
  etc. work; keep the cluster key `leiden`.
