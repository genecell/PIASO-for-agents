# Workflow — streaming / out-of-core analysis on a `.cytome`

The same pipeline as `end_to_end_scrnaseq.md`, run **from a file on disk**: peak memory is set by
`batch_size`, not by cell count; results (metrics, clusters, embeddings, labels, scores) are
written back onto the file and persist across sessions; an R collaborator opens the same file with
`read_cytome()`. Executed on the `e18_v3_nuclei` fixture converted to a cytome; the piaso.org
twins (*… end to end (cytome)*, *GDR at scale: 200,000 cells in 17 minutes*, *GDR and SVD on 1.5
million cells*) run the identical calls on atlases.

## Install

```bash
pip install piaso-tools      # cytome and cosg come with it
```

## When to convert (decision rule)

More than ~100k cells; several samples in one object; ATAC fragments; a tissue image with ROI
queries; results that should live on disk; an R handoff. Otherwise stay in AnnData. Either way the
analysis code is identical — only the first argument changes.

## Step 1 — Convert once (three entry points)

```python
import numpy as np, pandas as pd
import cytome, piaso, cosg

# a) Cell Ranger h5 -> cytome, no AnnData in between (RNA + ATAC for multiome)
ds = cytome.from_10x_h5("e18_v3_nuclei.h5", "e18.cytome", sample_name="E18")   # returns the OPEN dataset; keep using it
ds.n_cells, ds.n_genes, sorted(ds.modalities)                                    # (5973, 32285, ['RNA'])

# b) an .h5ad too large to load: stream it in
# ds = cytome.from_h5ad("big.h5ad", output="big.cytome", modality="RNA", backed=True, chunk_size=2048)

# c) an AnnData already in memory (keeps layers, obsm as RNA_umap etc.)
# ds = cytome.from_anndata(adata, modality="RNA", output="e18.cytome")

# several samples -> one file
# ds = cytome.merge(["s1.cytome", "s2.cytome"], output="merged.cytome", batch_key="sample_id")
```
**One writer per file**: use the returned Dataset; do not `cytome.open()` the same path while it
is open. `ds.flush()` after your own writes, `ds.close()` at the end, `cytome.open(path)` to resume.

## Step 2 — QC by streaming

```python
piaso.pp.calculateCellMetrics(ds, modality="RNA", prefix_vars={"mt": "mt-", "ribo": ["Rps", "Rpl"]})
cells = ds.cells.to_pandas()                        # the cells table is small; the MATRIX stays on disk
cells[["n_counts", "n_genes", "pct_counts_mt"]].describe()
piaso.pl.scatter(ds, x="n_counts", y="n_genes", color="pct_counts_mt", logx=True, logy=True, marginals=True)
piaso.pp.scrublet(ds, expected_doublet_rate=0.06, random_state=0)   # -> cells['scrublet_score', 'is_doublet']
```
Metrics are computed in cell batches and written into the `cells` table — on the file, not in a
variable.

## Step 3 — Normalize, embed, cluster (identical calls)

```python
piaso.tl.infog(ds, modality="RNA", n_top_genes=3000, save_layer=True)   # save_layer=True persists RNA_infog (disk); False recomputes on the fly
piaso.tl.runSVD(ds, modality="RNA", layer="infog", n_components=50, key_added="X_svd")
piaso.tl.neighbors(ds, use_rep="X_svd", n_neighbors=15)
piaso.tl.leiden(ds, resolution=1.0, key_added="leiden")                # -> cells['leiden']
piaso.tl.umap(ds, use_rep="X_svd")                                     # -> embedding 'X_umap'
piaso.pl.embedding(ds, basis="X_umap", color="leiden", legend_loc="both")
piaso.pl.plot_features_violin(ds, ["n_genes", "n_counts", "pct_counts_mt", "scrublet_score"], groupby="leiden")
```

## Step 4 — Markers: the one call that returns instead of writing

```python
markers = cosg.cosg(ds, groupby="leiden", modality="RNA", n_genes_user=25, layer="infog")   # dict; no key_added (no .uns)
order = list(markers["groups_order"])
names = pd.DataFrame(markers["names"], columns=order)
scores = pd.DataFrame(markers["scores"], columns=order)

qc = cells.groupby("leiden", observed=True).agg(median_scrublet=("scrublet_score", "median"),
        frac_high_scrublet=("scrublet_score", lambda v: float((v > 0.2).mean())))
qc = qc.loc[[c for c in order if c in qc.index]]                   # ALIGN to COSG's order before zipping
qc["top_cosg"] = [float(scores[c].iloc[0]) for c in qc.index]
qc["mt_in_top10"] = [sum(str(g).startswith("mt-") for g in names[c].head(10)) for c in qc.index]
drop = list(qc.index[(qc["median_scrublet"] > 0.15) | (qc["frac_high_scrublet"] > 0.5) | (qc["mt_in_top10"] >= 3)])

top3 = []
for c in order:
    top3 += [g for g in map(str, names[c][:3]) if g not in top3]
piaso.pl.dotplot(ds, top3, groupby="leiden", modality="RNA", cytome_layer="infog", standard_scale="var")   # reads only the plotted genes
```
Re-index per-cluster tables by `groups_order`: the groupby sorts labels as strings, COSG returns
its own order, and zipping them attaches each cluster's markers to another cluster's statistics —
plausibly.

## Step 5 — Filter writes a new file, then re-run

```python
cells = ds.cells.to_pandas()
hi_c = np.percentile(cells["n_counts"], 99)
keep = (~cells["leiden"].astype(str).isin(drop) & (cells["n_counts"] >= 500) & (cells["n_counts"] <= hi_c)
        & (cells["n_genes"] >= 250) & ~(cells["is_doublet"].astype(bool) | (cells["scrublet_score"] > 0.3)))
piaso.pp.filter_cells(ds, mask=np.asarray(keep), inplace=False, output="e18_clean.cytome", overwrite=True)
ds.close(); ds = cytome.open("e18_clean.cytome")

piaso.tl.infog(ds, modality="RNA", n_top_genes=3000, save_layer=True)
piaso.tl.runSVD(ds, modality="RNA", layer="infog", n_components=50, key_added="X_svd")
piaso.tl.neighbors(ds, use_rep="X_svd", n_neighbors=15); piaso.tl.leiden(ds, resolution=1.0, key_added="leiden"); piaso.tl.umap(ds, use_rep="X_svd")
```
A cytome subset is a copy: `inplace=False` + `output=` leaves the original intact;
`inplace=True, output=None` replaces the file atomically. `cell_mask=` on any plot/function selects
without writing.

## Step 6 — Annotation and scores, written back to the file

```python
markers_df, marker_sets = piaso.tl.getMarkers(study="AllenWholeMouseBrain_isocortex", as_dict=True)
piaso.tl.predictCellTypeByMarker(ds, marker_gene_set=marker_sets, modality="RNA", cytome_layer="infog",
                                 score_layer="infog", use_rep="X_svd", key_added="CellTypes")   # -> cells['CellTypes']
ds.cells.to_pandas()["CellTypes"].value_counts().head()
piaso.pl.embedding(ds, basis="X_umap", color="CellTypes"); piaso.pl.stackedBarplot(ds, groupby="CellTypes", splitby="leiden")

piaso.tl.score(ds, gene_list=["Gad1", "Gad2", "Slc32a1"], key_added="gaba", modality="RNA", layer="infog")   # -> cells['gaba']
piaso.tl.runGDR(ds, groupby="CellTypes", layer="infog", modality="RNA", key_added="X_gdr")
piaso.tl.neighbors(ds, use_rep="X_gdr", n_neighbors=15, key_added="gdr"); piaso.tl.umap(ds, use_rep="X_gdr", key_added="X_umap_gdr", neighbors_key="gdr")
ds.set_categories("CellTypes", order=sorted(ds.cells.to_pandas()["CellTypes"].astype(str).unique()))   # order + palette persist in the file
ds.close()
```
Nothing has to be exported for the annotation to persist: it is on the file the moment the call
returns, and `cytome.open("e18_clean.cytome")` in a new session (or `read_cytome()` in R) sees it.

## At atlas scale

The piaso.org *GDR at scale* run is the same six calls on a **path** (`piaso.tl.infog(path, ...)`,
`runGDR(path, batch_key="library_prep", groupby=None, ...)`) — 200,061 cells in 17 min, 1.5 M cells
with peak memory set by `batch_size`. Registry atlases: `piaso.data.load_dataset("allen_devvis_rna",
return_type="cytome")` (1.4 GB), `"humanlifespan_pfc_rna"` (25.7 GB). Knobs: `batch_size`,
`max_workers`, `max_score_batch_cache_bytes`.

## Handoff to R

```r
library(cytome)
so  <- read_cytome("e18_clean.cytome", as = "Seurat")          # clusters, CellTypes, embeddings, graphs come along
sce <- read_cytome("e18_clean.cytome", delayed = TRUE)         # out-of-core for scran / scuttle
```
And back: `write_cytome(so, "from_r.cytome")` → any Python call above. See `components/cytome.md`.

## Back to AnnData

`adata = ds.to_anndata(modality="RNA")` restores `obsm['X_umap']`, `obsm['X_gdr']`,
`layers['infog']` and the `cells` columns as `obs` — for anything that still wants an in-memory
object.
