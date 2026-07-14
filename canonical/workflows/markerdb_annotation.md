# Workflow — cell-type inference from PIASOmarkerDB

Infer cell types for a gene list by querying **PIASOmarkerDB**, a curated marker database, over its
live REST API. Use this when you have a set of genes (e.g. cluster markers) and want to know which
cell types they point to, drawing on published atlases rather than your own reference. Both blocks
are executed and passing.

## Install

```bash
pip install piaso-tools "matplotlib<3.9"
```

**Internet required.** PIASOmarkerDB is a **remote REST API client** (base
`https://piaso.org/piasomarkerdb/api/v1/`) — it is not bundled data, so these blocks need network
egress to `piaso.org`. This workflow is **Python-only**; there is no R client.

## Step 1 — Query the database

Explore what the DB contains and pull markers by filter. Returns pandas DataFrames (or lists for the
`list_*` meta-queries). No AnnData involved.

```python
import piaso
studies = piaso.tl.queryPIASOmarkerDB(list_studies=True)      # list[str], 36 studies
df = piaso.tl.queryPIASOmarkerDB(cell_type="Microglia", limit=5)
# columns: cell_type, condition, gene, species, specificity_score, study_publication, tissue
df2 = piaso.tl.getMarkers(gene="AIF1", limit=3)               # getMarkers is an alias of queryPIASOmarkerDB
```
**Out:** DataFrames of marker records; `studies` is the list of available study keys. Filters
accepted: `gene`, `cell_type`, `study`, `species`, `tissue`, `condition`, `min_score`, `max_score`,
`limit`.

## Step 2 — Infer cell types for a gene list

Pass a gene list to `analyzeMarkers`; it queries the DB for each gene, groups hits by
(cell_type, study, species, tissue, condition), and ranks contexts by matched-gene count then
average specificity. A plain list returns a ranked DataFrame.

```python
res = piaso.tl.analyzeMarkers(["AIF1", "P2RY12", "CX3CR1", "CSF1R"], n_top_genes=50)  # infer cell type
# DataFrame cols: cell_type, study_publication, species, tissue, condition,
#                 matched_gene_count, matched_genes, avg_specificity
```
**Out:** ranked DataFrame of candidate cell types for the gene list (top row = best match).

## Notes

- `analyzeMarkers` also accepts a COSG-style DataFrame (columns = clusters) or a
  `{cluster: [genes]}` dict; those inputs return a `(results_dict, top_hits)` **tuple**, where
  `top_hits[cluster]` is a plain **cell-type-name string** (or `"Unassigned"`) — not a DataFrame
  or nested dict. Map it straight onto clusters, e.g.
  `adata.obs["cell_type"] = adata.obs["Leiden"].map(top_hits)`.
- This composes with `marker_based_annotation.md`: `predictCellTypeByMarker` transfers the
  *marker-set keys* onto cells (cluster IDs, if that's what you keyed the marker set by), so to
  attach real biological names, feed the same COSG marker set to `analyzeMarkers` and map its
  `top_hits` back onto the clusters.
- Narrow noisy results with `species` / `tissue` / `studies` (validated against `list_studies`;
  unknown study names raise a `ValidationError`), or `exclude_studies` / `exclude_cell_types`.
- This complements `marker_based_annotation.md`: that workflow transfers **your** cluster labels via
  scoring; this one names cell types from a **curated public** database.
