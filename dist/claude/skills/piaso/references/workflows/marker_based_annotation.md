# Workflow — marker-based cell-type annotation (and reference projection)

Annotate cell types by scoring each cell against per-type marker sets with
`piaso.tl.predictCellTypeByMarker`, smoothing over an embedding. Three ways to get the marker
sets and the embedding, in order of preference; all blocks are executed on the `e18_v3_nuclei`
fixture (see `data.md`). Mirrors the executed piaso.org tutorials *Marker-based cell type
prediction* (95.4 % held-out accuracy on the cortex reference) and *projectGDR*.

## Install

```bash
pip install piaso-tools
```

## Prerequisites

`adata` from `end_to_end_scrnaseq.md` through Step 7: `layers["infog"]`, `obs["leiden"]`,
`obsm["X_svd"]` (and `uns["cosg"]` if you ran COSG). No `X_gdr` is required — build it below.

## Decision rule — which route

| You have | Route | Function |
|---|---|---|
| A matching curated study (tissue + species) | **A** — PIASOmarkerDB marker sets | `getMarkers(study=, as_dict=True)` → `predictCellTypeByMarker` |
| An annotated reference dataset of the same tissue | **B** — COSG markers from the reference; optionally its GDR space | `cosg.cosg(ref, groupby=)` → dict → `predictCellTypeByMarker`; or `runGDR(ref, save_reference=True)` + `projectGDR(query, ref)` |
| Only your clusters' markers | **C** — name them | `analyzeMarkers(dict)` → `top_hits` (see `markerdb_annotation.md`) |
| A reference and you want a joint embedding | **D** | `predictCellTypeByGDR(query, ref, ...)` (≥ 1.2.3) |

## Route A — curated marker sets (PIASOmarkerDB)

```python
import pandas as pd
import piaso, cosg
piaso.tl.getMarkers(list_studies=True)                     # 36 studies; pick the tissue/species match
markers_df, marker_db = piaso.tl.getMarkers(study="AllenWholeMouseBrain_isocortex", as_dict=True)   # TUPLE — unpack both
len(marker_db)                                             # 26 cell types -> gene lists

piaso.tl.predictCellTypeByMarker(adata, marker_gene_set=marker_db, score_layer="infog",
                                 use_rep="X_svd", key_added="CellTypes_db")
# writes obs['CellTypes_db'] + CellTypes_db_raw / _smoothed / _score / _smoothed_confidence / _confidence_smoothed
# and obsm['CellTypes_db_score'] (cells x types)
```
`use_rep` is the embedding used for kNN smoothing (`k_nearest_neighbors=7`); its default is `X_gdr`,
so pass `X_svd` when that is what you have. `score_method="piaso"` (PIASOscore) is the default.

## Route B — an annotated reference

```python
ref = piaso.data.load_dataset("adult_cortex_multiome_rna")          # 17,412 cells, 20 curated types, 5 protocols; counts in layers['raw']
piaso.tl.infog(ref, layer="raw", n_top_genes=3000)                  # infog REFUSES the scaled .X and names the layer
cosg.cosg(ref, groupby="CellTypes", key_added="cosg", n_genes_user=30, mu=10, layer="infog")   # mu=10: specific markers (Pvalb gains Tac1, Syt2)
names = pd.DataFrame(ref.uns["cosg"]["names"])
marker_ref = {ct: list(names[ct]) for ct in names.columns}

piaso.tl.predictCellTypeByMarker(adata, marker_gene_set=marker_ref, score_layer="infog",
                                 use_rep="X_svd", key_added="CellTypes_ref")
```

### B′ — the reference's own GDR space, projected

When the query should live in the reference's coordinates (atlas + incoming samples, gates that
must stay valid), freeze the reference's GDR and project:

```python
piaso.tl.runGDR(ref, batch_key="Sample", groupby=None, layer="infog", infog_layer="raw",
                score_layer="infog", n_gene=30, key_added="X_gdr", save_reference=True)   # uns['gdr_reference']: marker sets + norms
piaso.tl.projectGDR(adata, reference=ref, key_added="X_gdr_proj")   # query gets the REFERENCE's width; streams on a .cytome query
piaso.tl.predictCellTypeByMarker(adata, marker_gene_set=marker_ref, score_layer="infog",
                                 use_rep="X_gdr_proj", key_added="CellTypes_ref")
```
`mode="reference"` (default) keeps coordinates comparable across queries; `mode="self"` treats the
query as a new batch (falls back below 500 cells). `novelty_k` / `novelty_quantile` flag cells far
from any reference cell — a type the reference lacks gets coordinates but not meaning. Executed on
the fixture: `projectGDR` on a 30 % hold-out reproduces the reference's 22-dimensional space.

## Route D — joint embedding (`predictCellTypeByGDR`, ≥ 1.2.3)

```python
piaso.tl.predictCellTypeByGDR(query, ref, layer="infog", layer_reference="infog",
                              reference_groupby="CellTypes", query_groupby="leiden", key_added="CellTypes_gdr")
```
Concatenates reference + query, runs GDR jointly, Harmony-integrates, trains an SVM on the
reference, predicts the query (`obs["CellTypes_gdr"]`). Query and reference must have **distinct
`obs_names`** (disjoint cells), and the `reference_groupby` column must exist on the query too
(`query.obs["CellTypes"] = "unknown"` as a placeholder). Prefer B′ when the reference's coordinates must not move.

## Check the annotation, then use it

```python
# do the two sources agree?
from sklearn.metrics import adjusted_rand_score
adjusted_rand_score(adata.obs["CellTypes_db"], adata.obs["CellTypes_ref"])     # 0.88 on the mouse-brain tutorial; broad classes agree, fine excitatory subtypes do not

# markers that drove the call, grouped by the label it produced — a clean diagonal is evidence
top, seen = [], set()
for ct in adata.obs["CellTypes_db"].astype("category").cat.categories:
    picked = [g for g in marker_db[ct] if g in adata.var_names and g not in seen][:3]
    top += picked; seen.update(picked)
piaso.pl.dotplot(adata, top, groupby="CellTypes_db", standard_scale="var")
piaso.pl.plotConfusionMatrix(adata, groupby_query="CellTypes_db", groupby_reference="CellTypes_ref")
piaso.pl.sankey(adata, left="CellTypes_ref", right="CellTypes_db")

# per-type QC after annotation: a type that sits apart on a technical metric was called on that metric
gm = piaso.pp.calculateGroupMetrics(adata, groupby="CellTypes_db")
piaso.pl.plotGroupMetrics(gm, data=adata, groupby="CellTypes_db")

# rebuild the embedding from the labels
piaso.tl.runGDR(adata, groupby="CellTypes_db", layer="infog", key_added="X_gdr")
```

## Notes

- `marker_gene_set` may be any `{label: [genes]}` dict, a DataFrame (one column per type) or a
  list — you are not limited to COSG or PIASOmarkerDB. Labels are the dict's **keys**.
- `getMarkers(..., as_dict=True)` **returns a tuple**; assigning it to one name hands the tuple to
  `predictCellTypeByMarker`. Cell-type names follow the study's vocabulary exactly.
- Held-out truth on the cortex reference: 95.4 % with `mu=10` markers on a GDR embedding vs
  91.5 % with `mu=1` markers on SVD — the embedding and marker specificity are worth ~4 points
  together; the residual errors are adjacent cortical layers.
- On a `.cytome`: same calls with `modality="RNA"`, `cytome_layer="infog"`; labels land in the
  `cells` table (`streaming_large_data.md`).
