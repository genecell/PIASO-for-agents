# Workflow — cell-type inference from PIASOmarkerDB

Query **PIASOmarkerDB**, a curated marker database (36 studies, human + mouse; brain, blood, bone
marrow, breast, thymus, whole-body atlases), over its live REST API: pull marker sets by study /
gene / cell type, or go the other way — from a gene list (a cluster's COSG markers, a DE hit list)
to the cell types it matches. All blocks are executed. Mirrors the executed piaso.org tutorial
*PIASOmarkerDB API client*.

## Install

```bash
pip install piaso-tools
```

**Internet required.** The client issues HTTP calls to `https://piaso.org/piasomarkerdb/api/v1/`
(cache under `~/.piaso/markers`); nothing is bundled. Python-only — from R or another agent use
the `piaso-mcp` server (`query_marker_db`, `list_studies`) or the REST endpoints
`/markers?gene=|cell_type=|study=|species=|tissue=|limit=`, `/studies`, `/genes` directly.

## Step 1 — What is in there, and one study

```python
import pandas as pd
import piaso
studies = piaso.tl.getMarkers(list_studies=True)                  # list[str], 36 studies (getMarkers == queryPIASOmarkerDB)
df = piaso.tl.getMarkers(study="AllenWholeMouseBrain_isocortex")  # (1300, 7): cell_type, condition, gene, species, specificity_score, study_publication, tissue
markers_df, marker_sets = piaso.tl.getMarkers(study="AllenWholeMouseBrain_isocortex", as_dict=True)   # TUPLE: table + {cell_type: [genes]} (26 types)
```
`specificity_score` is a COSG score (higher = more exclusive to that type). **`as_dict=True`
returns both** — unpack two names; the dict is what `predictCellTypeByMarker` consumes.

## Step 2 — Ask the other way round

```python
piaso.tl.getMarkers(gene="Sst")                                    # where is Sst a marker, across studies (interneurons AND gut enteroendocrine cells)
piaso.tl.getMarkers(cell_type="Pvalb Gaba")                        # names must match the study vocabulary exactly
piaso.tl.getMarkers(species="Mouse", limit=5)
piaso.tl.getMarkers(study="AllenWholeMouseBrain_isocortex", min_score=5.0)
piaso.tl.getMarkers(study="AllenWholeMouseBrain_isocortex", list_cell_types=True)[:6]   # ['001 CLA-EPd-CTX Car3 Glut', '004 L6 IT CTX Glut', ...]
len(piaso.tl.getMarkers(study="AllenWholeMouseBrain_isocortex", list_genes=True))      # 23657
```
Filters: `gene`, `cell_type`, `study`, `species`, `tissue`, `condition`, `min_score`, `max_score`,
`limit`. `"Astro-TE"` returns 0 rows because the stored name is `"319 Astro-TE NN"` — list first.

## Step 3 — From a gene list to a cell type

```python
res = piaso.tl.analyzeMarkers(["Sst", "Pvalb", "Vip", "Lamp5", "Gad1", "Gad2"])
# ranked DataFrame: cell_type, study_publication, species, tissue, condition, matched_gene_count, matched_genes, avg_specificity
# every top hit is an inhibitory neuron type across four studies and two species
```
Ranking is by matched-gene count, then average specificity. Narrow noisy results with
`species=`, `tissue=`, `studies=` (validated against `list_studies`), `exclude_studies=`,
`exclude_cell_types=`.

## Step 4 — Name your clusters

Feed the per-cluster COSG output (dict or the COSG DataFrame); the return becomes a **tuple**
`(results_dict, top_hits)` with `top_hits[cluster]` a cell-type **string** (`"Unassigned"` if nothing
matched):

```python
import cosg
cosg.cosg(adata, groupby="leiden", key_added="cosg", n_genes_user=25, layer="infog")
names = pd.DataFrame(adata.uns["cosg"]["names"])
results, top_hits = piaso.tl.analyzeMarkers({c: list(names[c]) for c in names.columns},
                                            n_top_genes=25, species="Mouse")
adata.obs["celltype_db"] = adata.obs["leiden"].map(top_hits)      # e.g. '0' -> '038 DG-PIR Ex IMN', '2' -> 'UL CPN'
results["0"].head()                                                # the ranked candidates behind each call
```
This is the fastest sanity check on an unnamed cluster. It **names** clusters; to label
**cells** by score (and smooth over the embedding) use `predictCellTypeByMarker` with the study's
`marker_sets` (`marker_based_annotation.md`). The two compose: predict cells from a study's sets,
then confirm cluster identities with `analyzeMarkers` on your own COSG markers.

## Offline alternative

`piaso.data.load_dataset("piaso_markerdb_allen_immune")` returns a 115 KB static slice
(Allen Human Immune Health Atlas L2) as a DataFrame — for tests and air-gapped work; it is a
published export, not the API.
