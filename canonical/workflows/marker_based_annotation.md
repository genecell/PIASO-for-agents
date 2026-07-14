# Workflow — marker-based cell-type annotation

Annotate cell types by scoring each cell against per-cluster marker sets. Builds on a clustered
AnnData: derive markers with COSG, then predict labels with PIASO's `predictCellTypeByMarker`,
smoothing over the marker-guided embedding. Both blocks are executed and passing.

## Install

```bash
pip install piaso-tools cosg "matplotlib<3.9"
```

## Prerequisites (what must already exist)

This workflow **continues from `end_to_end_scrnaseq.md`**. Before starting, `adata` must have:
`adata.layers["infog"]` (the scoring layer), `adata.obs["Leiden"]` (clusters), and
`adata.obsm["X_gdr"]` (used to smooth predictions). Run the end-to-end workflow through Step 6 first
if these are absent.

## Step 1 — COSG markers → marker set dict

Runs COSG on the clusters and reshapes its recarray output into a `{cluster: [genes]}` dict, which
is the input `predictCellTypeByMarker` expects. (If you already ran Step 5 of the end-to-end
workflow, `adata.uns["cosg"]` exists and you can skip straight to building `marker_set`.)

```python
import cosg
adata.X = adata.layers["infog"]          # COSG expects normalized/log values in .X (or a layer)
cosg.cosg(adata, groupby="Leiden", key_added="cosg", n_genes_user=30, mu=1.0)

names = adata.uns["cosg"]["names"]
marker_set = {cl: [names[cl][i] for i in range(len(names))] for cl in names.dtype.names}
```
**Out:** `marker_set` — `{cluster: [top genes]}`.

## Step 2 — Predict cell types from the marker set

Scores every cell against each cluster's marker set (over the `infog` layer), assigns the argmax
label, then smooths the prediction across the `X_gdr` neighborhood. Writes the prediction plus a
score matrix and confidence.

```python
import piaso
piaso.tl.predictCellTypeByMarker(adata, marker_gene_set=marker_set, score_layer="infog",
                                 use_rep="X_gdr", key_added="pred", smooth_prediction=True,
                                 return_confidence=True)
# writes adata.obs['pred'] (+ pred_raw/pred_smoothed/pred_score/pred_*_confidence) and adata.obsm['pred_score']
```
**Out:** `adata.obs["pred"]` (final label) + `adata.obsm["pred_score"]` (full score matrix). The
`.obs` also gets `pred_score`, `pred_raw`, `pred_smoothed`, `pred_smoothed_confidence`, and
`pred_confidence_smoothed` (verified against the current build — the exact set written).

> **The `pred` labels mirror the keys of `marker_gene_set`.** If you keyed the marker set by
> cluster ID (as above), `pred` is a (smoothed) cluster ID, not a biological name. To attach real
> cell-type names, feed the same COSG marker set to `analyzeMarkers` (see
> `markerdb_annotation.md`) and map its `top_hits` onto the clusters.

## Notes

- `marker_gene_set` can also be a curated set from PIASOmarkerDB or any hand-authored
  `{cell_type: [genes]}` dict — you are not limited to COSG output.
- `score_layer="infog"` and `use_rep="X_gdr"` are defaults that **must exist** on the object;
  smoothing is skipped or errors if `X_gdr` is missing.
- To infer cell-type *names* (rather than transfer cluster labels) from a marker gene list against a
  curated database, use `markerdb_annotation.md` instead.
