# PIASO — component reference

PIASO is a Python single-cell omics toolkit (Gord Fishell Lab, HMS / Broad). It bundles
INFOG normalization, SVD and marker-gene-guided dimensionality reduction (GDR), gene-set
enrichment scoring (Rust-accelerated), marker-based cell-type prediction, marker-guided
batch integration, single-cell ligand–receptor inference (SCALAR), and a client for the
remote PIASOmarkerDB marker database. This file is self-sufficient: it assumes nothing is
already installed or imported.

## Install

```bash
pip install piaso-tools "matplotlib<3.9"
```

The `matplotlib<3.9` pin is **mandatory**. PIASO 1.1.0 imports `piaso.plotting.color` at
package import time, and that module calls the `matplotlib.cm.get_cmap(...)` API that was
**removed in matplotlib 3.9**. Under matplotlib ≥ 3.9, `import piaso` fails outright.
PIASO's own `pyproject.toml` only requires `matplotlib>=3.5.2` (no upper cap), so the cap
is effectively undeclared — pin it yourself in every install line until fixed upstream.

The PyPI distribution is named `piaso-tools`; the import name is `piaso`. Last tested
against version 1.1.0.

## Import convention

Every public symbol lives under a submodule — there is **no** top-level `piaso.X`
re-export. Use the short aliases consistently:

```python
import piaso
piaso.tl   # tools          (also available as piaso.tools)
piaso.pp   # preprocessing  (also available as piaso.preprocessing)
piaso.pl   # plotting       (also available as piaso.plotting)
```

The short and long forms are the **same runtime objects** (`piaso.tools is piaso.tl`
returns `True`), so `piaso.tl.runGDR` and `piaso.tools.runGDR` are identical. This file
uses the short forms throughout.

## Dependency: COSG

`piaso-tools` **hard-depends on `cosg`** (auto-installed with it). PIASO calls
`cosg.cosg(...)` internally inside GDR and stitchSpace, so those functions fail if COSG is
missing — but COSG is **not re-exported** under the `piaso.*` namespace. There is no
`piaso.cosg`; to call COSG directly you still `import cosg`. See `components/cosg.md`.

## Citation

> Wu, S.J., Dai, M. et al. Pyramidal neurons proportionately alter cortical interneuron
> subtypes. Nature (2026). DOI: 10.1038/s41586-025-09996-8

---

## Shared setup for the examples

The code blocks below all build on a loaded, filtered AnnData with raw UMI counts kept in a
`counts` layer:

```python
import scanpy as sc
adata = sc.read_10x_h5("e18_v3_nuclei.h5")   # any 10x h5; ~5k cells
adata.var_names_make_unique()
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.layers["counts"] = adata.X.copy()      # keep raw UMIs
```

---

## INFOG normalization — `piaso.tl.infog`

**What it computes:** an information-content-based normalization of **raw UMI counts** (not a
scoring or DR method). It depth-normalizes counts to the median library size, scales by an
information factor, takes an elementwise square root, optionally trims extreme values, then
selects the top-`n_top_genes` highly variable genes by variance of the normalized matrix.

**Reads:** raw UMI counts from `adata.X` or `adata.layers[layer]` (negatives are rejected).
**Writes:** normalized matrix to `adata.layers[key_added]` (default `infog`) — or to
`adata.X` when `inplace=True`; plus `adata.var['<key_added>_var']` and
`adata.var['highly_variable']` (bool).

INFOG must be given **raw counts**. Its output layer name defaults to `infog`, which is also
the default input layer of `piaso.tl.score` (below).

```python
import piaso
piaso.tl.infog(adata, layer=None, n_top_genes=2000, key_added="infog")
# writes adata.layers['infog'] (normalized) + adata.var['highly_variable']
```

---

## SVD embedding — `piaso.tl.runSVDLazy` / `piaso.tl.runSVD`

**What they compute:** a low-dimensional cell embedding via sklearn `TruncatedSVD`
(randomized). `runSVD` is the core call (assumes an `highly_variable` column already
exists). `runSVDLazy` is the workhorse wrapper: it does HVG selection + SVD in one call, and
when `layer='infog'` it runs INFOG normalization first (from a raw-counts layer) and does
SVD on that.

**`runSVDLazy` reads:** `adata.X` / `adata.layers[layer]`; for infog mode, raw counts from
`adata.layers[infog_layer]`. **Writes:** `adata.obsm[key_added]` (default `X_svd`, shape
`n_obs × n_components`); in infog mode also `adata.layers['infog']`,
`adata.var['infog_var']`, `adata.var['highly_variable']`.

Note the default `random_state` differs between the two: `runSVDLazy` uses 1927, `runSVD`
uses 10. `runSVDLazy` is what GDR, `leiden_local`, and the parallel wrappers use internally.

```python
piaso.tl.runSVDLazy(adata, layer="infog", infog_layer="counts",
                    n_components=50, n_top_genes=2000, key_added="X_svd", random_state=1927)
# writes adata.obsm['X_svd'] (n_obs x 50)
```

The embedding then feeds a standard neighbors + Leiden clustering step (scanpy; needs
`igraph` + `leidenalg`):

```python
sc.pp.neighbors(adata, use_rep="X_svd", n_neighbors=15)
sc.tl.leiden(adata, resolution=1.0, key_added="Leiden", flavor="igraph",
             n_iterations=2, directed=False)
```

---

## Marker-gene-guided DR (GDR) — `piaso.tl.runGDR`

**What it computes:** GDR is **marker Gene-guided Dimensionality Reduction**, not a generic
matrix factorization. It (1) takes cluster labels (`groupby`, or de-novo clustering);
(2) runs **COSG** to get the top-`n_gene` marker genes per cluster; (3) scores every cell
against each cluster's marker set (scanpy `score_genes` by default, or PIASO `score` when
`scoring_method='piaso'`); (4) double L2-normalizes the resulting cell × cluster score
matrix. **That marker-score matrix is the embedding.**

**Key consequence:** the embedding width `X_gdr` is the **number of clusters**, not
`n_svd_dims`. With a `batch_key`, markers are found per batch and every cell is scored
against all batches' marker sets, then horizontally stacked → a batch-integrated embedding.

**Reads:** `adata.X` / `adata.layers[layer]` (log/normalized), `adata.obs[groupby]` (and
`[batch_key]`), optional `score_layer` / `infog_layer`. Requires `cosg` installed.
**Writes:** `adata.obsm[key_added or 'X_gdr']` (`n_obs × total_n_clusters`) and
`adata.uns['gdr']`. Its default `scoring_method` resolves to `'scanpy'` (note
`runGDRParallel` instead defaults to `'piaso'`).

The GDR marker step needs COSG-computed markers on normalized values. Run COSG first
(COSG expects normalized/log values in `.X` or a layer):

```python
import cosg
adata.X = adata.layers["infog"]          # COSG expects normalized/log values in .X (or a layer)
cosg.cosg(adata, groupby="Leiden", key_added="cosg", n_genes_user=30, mu=1.0)
# writes adata.uns['cosg'] = {'names': recarray, 'scores': recarray}, one field per cluster
```

```python
piaso.tl.runGDR(adata, groupby="Leiden", n_gene=30, mu=1.0, layer="infog",
                score_layer="infog", scoring_method="scanpy", key_added="X_gdr")
# writes adata.obsm['X_gdr'] (n_obs x n_clusters)
```

The output `X_gdr` is used as `use_rep` for neighbors/UMAP/Leiden and by
`predictCellTypeByMarker` / `predictCellTypeByGDR`.

---

## Gene-set enrichment scoring — `piaso.tl.score`

**Source-only / UNDOCUMENTED** — this function is not on the PIASO website (a genuine
selling point to surface). It is **Rust-accelerated**: the multi-gene-set path lazily calls
the bundled Rust extension `piaso._piaso_score` (`score_complete` / `fused_matmul_reduce`,
which release the GIL), with a pure-Python fallback when the extension is absent.

**What it computes:** gene-set enrichment scoring with expression-matched control-gene
background subtraction — PIASO's own optimized algorithm and implementation. For each gene
set it builds control sets by KNN in (mean, variance) space, computes the weighted query
score minus the mean control score, plus empirical/Monte-Carlo p-values and BH-FDR. Behavior
depends on the input type:

- **Single set** (`list[str]`): writes to `adata`, returns `None`.
- **Multiple sets** (`dict` / `DataFrame` / `list[list]`): returns
  `(score_matrix, gene_set_names, pval_matrix)` and uses the Rust backend.

**Reads:** `adata.layers[layer]` (default `'infog'`) or `adata.X`; expects INFOG-normalized
values by default. **Writes (single-set only):** `adata.obs[key_added or 'INFOG_score']`
plus `adata.uns[key_added or 'INFOG_score']` (a DataFrame of score / query / ctrl / pvals).
Genes not in `var_names` are silently dropped; the default `layer='infog'` errors if that
layer is absent.

```python
genes = [g for g in ["Gad1","Gad2","Slc17a7"] if g in adata.var_names]
piaso.tl.score(adata, gene_list=genes, layer="infog", key_added="myscore")
# writes adata.obs['myscore'] + adata.uns['myscore'] (score/query/ctrl/pvals)
```

---

## Cell-type prediction — `piaso.tl.predictCellTypeByMarker` / `piaso.tl.predictCellTypeByGDR`

Two complementary annotation routes.

### `predictCellTypeByMarker` — marker-set based

**What it computes:** scores each cell against every cell type's marker set (via the parallel
scorer), predicts the label as the argmax score, then optionally smooths predictions over a
kNN graph in `use_rep`. Input marker sets can come from COSG output or from PIASOmarkerDB.

**Reads:** `adata.X` / `adata.layers[score_layer]` (default `infog`), the `marker_gene_set`
(list/dict/DataFrame), and `adata.obsm[use_rep]` (default `X_gdr`) for smoothing — both
must exist. **Writes:** `adata.obs[key_added]` (final label), `adata.obsm[key_added+'_score']`
(full score matrix), plus confidence and `_smoothed` / `_raw` variants.

```python
names = adata.uns["cosg"]["names"]
marker_set = {cl: [names[cl][i] for i in range(len(names))] for cl in names.dtype.names}
piaso.tl.predictCellTypeByMarker(adata, marker_gene_set=marker_set, score_layer="infog",
                                 use_rep="X_gdr", key_added="pred", smooth_prediction=True,
                                 return_confidence=True)
# writes adata.obs['pred'] (+ pred_raw/pred_smoothed/pred_score/pred_*_confidence) and adata.obsm['pred_score']
```

### `predictCellTypeByGDR` — reference-based label transfer

**What it computes:** concatenates a reference and query AnnData, runs GDR on the combined
object, applies Harmony integration on `X_gdr`, then trains an RBF SVM on the reference cells
and predicts query labels. **Reads:** `adata.layers[layer]` / `adata_ref.layers[layer_reference]`
(overwrites `.X`), `adata_ref.obs[reference_groupby]`, `adata.obs[query_groupby]`.
**Writes:** `adata.obs[key_added or 'CellTypes_gdr']`. Requires Harmony
(`scanpy.external`) and overwrites `.X` with the chosen layer.

> No verified code block exists for `predictCellTypeByGDR` — signature/behavior above is
> from source, not from an executed run.

---

## Local sub-clustering — `piaso.tl.leiden_local`

**What it computes:** sub-clusters existing groups **locally**. For each selected group (or
all groups jointly), it subsets those cells, recomputes DR (`X_pca` via runSVDLazy, or a GDR
/ harmony variant), rebuilds neighbors, and runs Leiden, relabelling as `{group}-{local}`.

**Reads:** `adata.obs[groupby]`, `adata.X` (and `batch_key` for harmony variants).
**Writes:** `adata.obs[key_added]` (default `Leiden_local`, categorical). `dr_method` is
restricted to a fixed set (raises otherwise); harmony variants require a `batch_key`.

```python
piaso.tl.leiden_local(adata, groupby="Leiden", key_added="Leiden_local",
                      resolution=0.25, dr_method="X_pca")
# writes adata.obs['Leiden_local'] (categorical, '{group}-{local}' labels)
```

> No verified code block exists for `leiden_local` — the snippet above follows the source
> signature but was not part of the executed test set.

---

## Marker-guided batch integration — `piaso.tl.stitchSpace`

**IMPORTANT: `stitchSpace` is NOT spatial.** Despite the name, it reads no coordinates. It is
a **marker-gene-guided batch correction of an embedding** ("Space" = embedding space).

**What it computes:** Stage 1 builds a BBKNN graph across batches on `use_rep`,
Leiden-clusters within each batch, runs COSG per batch for markers, computes pairwise Jaccard
marker overlap between inter-batch clusters, and prunes BBKNN edges between incompatible
clusters. Stage 2 computes each cell's pruned-neighbor centroid and applies a single
correction step to the embedding.

**Reads:** `adata.obsm[use_rep]` (rejected if it contains NaN/Inf), `adata.obs[batch_key]`,
`adata.X` / `layers[filter_cosg_layer]` for COSG. **Writes:** `adata.obsm[key_added]`
(corrected embedding, same dim as `use_rep`), plus `adata.uns[...]` marker/param entries and
pruned-graph `adata.obsp[...]`. Requires `bbknn` and `cosg`.

> No verified code block exists for `stitchSpace`. It also carries the annoy segfault gotcha
> (see `gotchas.md`): `bbknn` with `annoy >= 1.17` can segfault — pin `annoy==1.16.3`.

---

## Spatial coordinate rotation — `piaso.pp.rotateSpatialCoordinates`

**The genuinely spatial helper in PIASO.** **What it computes:** a 2D rotation of spatial
coordinates about their centroid (extra dimensions pass through unchanged).

**Reads:** `adata.obsm[spatial_key]` (default `X_spatial`; KeyError if absent, ValueError if
< 2D). **Writes:** rotated coordinates back to that key (and to `backup_spatial_key` if
given); `inplace=False` returns a copy. Spatial transcriptomics only.

```python
piaso.pp.rotateSpatialCoordinates(adata, angle_degrees=90, spatial_key="X_spatial",
                                  clockwise=False, inplace=True)
# rotates adata.obsm['X_spatial'] in place
```

> No verified code block exists for `rotateSpatialCoordinates` — the snippet follows the
> source signature and needs a spatial fixture to test.

---

# SCALAR — single-cell ligand–receptor (`piaso.tl.runSCALAR`)

**SCALAR is a PIASO function, not a separate package.** It is `piaso.tl.runSCALAR`, shipped
inside `piaso-tools`. Do not look for a `scalar` PyPI package.

**What it computes:** cell-type-resolved ligand–receptor interaction inference for
**dissociated single-cell** (non-spatial) data. The interaction score is
`(ligand specificity in the sender cell type) × (receptor specificity in the receiver cell
type)`, read from a user-supplied specificity matrix. Significance comes from a
gene-expression-matched permutation null (control ligand/receptor genes sampled from a
per-gene KNN in mean/variance space), giving an empirical p-value and BH-FDR per
sender–receiver pair.

**User must supply the reference data — there is no bundled DB:**
- `specificity_matrix`: a genes × cell-types DataFrame. The tested recipe is a z-scored
  per-cluster mean (see the block below). COSG scores can seed it too, but `cosg.cosg` stores
  only the top-N genes per cluster in `adata.uns['cosg']`, so you'd have to assemble a full
  genes × cell-types matrix first — the z-scored mean is the simpler, all-genes option.
- `lr_pairs`: a DataFrame with `ligand` / `receptor` columns (optional pathway annotation
  column); genes should be present in `adata`. No bundled list, but if `laris` is installed you
  can reuse its CellChatDB: `laris.datasets.lrDatabase(species="mouse")[["ligand","receptor"]]`.

**Reads:** `adata.var_names`, `adata.X` / `adata.layers[layer]` (for the background
mean/variance KNN), and `specificity_matrix.index`. **Writes:** nothing to `adata` — it
**returns a DataFrame** with columns `ligand, receptor, sender, receiver, interaction_score,
p_value, p_value_fdr, nlog10_p_value_fdr`.

```python
import numpy as np, pandas as pd
# specificity_matrix: genes x cell types (z-scored per-cluster mean; tested recipe)
ct = adata.obs["Leiden"].astype(str)
cell_types = sorted(ct.unique(), key=int)
X = adata.layers["infog"]
means = np.vstack([np.asarray(X[(ct==c).values].mean(0)).ravel() for c in cell_types]).T
spec = pd.DataFrame(means, index=adata.var_names, columns=cell_types)
spec = spec.sub(spec.mean(1), axis=0).div(spec.std(1).replace(0,1), axis=0)
# lr_pairs: DataFrame with ligand/receptor columns, genes present in adata
lr = pd.DataFrame([{"ligand":l,"receptor":r} for l,r in
                   [("Nrxn1","Nlgn1"),("Nrxn3","Nlgn1"),("Efna5","Epha4")]
                   if l in adata.var_names and r in adata.var_names])
res = piaso.tl.runSCALAR(adata, specificity_matrix=spec, lr_pairs=lr,
                         n_permutations=200, random_seed=42)
# DataFrame: ligand, receptor, sender, receiver, interaction_score, p_value, p_value_fdr, nlog10_p_value_fdr
```

**Companion plotters** (consume the `runSCALAR` DataFrame + specificity matrix):
`piaso.pl.plotLigandReceptorInteraction` (barplot + specificity heatmap) and
`piaso.pl.plotLigandReceptorLollipop` (lollipop of top interactions).

**Decision rule — SCALAR vs LARIS.** For ligand–receptor / cell–cell communication:
- **Spatial data** (coordinates in `.obsm['spatial']` / `['X_spatial']`; Visium / MERFISH /
  Xenium) → use **LARIS** (`laris.tl.runLARIS`; bundles CellChatDB). See
  `components/laris.md`.
- **Dissociated single-cell** (no coordinates) → use **SCALAR** (`piaso.tl.runSCALAR`;
  user supplies the LR-pair list + specificity matrix).

---

# PIASOmarkerDB — remote marker database client

**PIASOmarkerDB is a remote REST API client, not bundled data.** The functions issue HTTP
queries to `https://piaso.org/piasomarkerdb/api/v1/` and therefore **require internet
access** (and the `requests` package; results cache under `~/.piaso/markers`). It is
**Python-only**.

Functions (all under `piaso.tl`):
- `queryPIASOmarkerDB(...)` — query the DB by gene / cell_type / study / species / tissue /
  condition / score range; also `list_studies` / `list_cell_types` / `list_genes`
  meta-queries. Returns a DataFrame (or a list for the `list_*` calls).
- `getMarkers(...)` — a thin `@wraps` **alias** of `queryPIASOmarkerDB` (identical behavior).
- `analyzeMarkers(genes, ...)` — infers likely cell types for a gene list / COSG DataFrame /
  dict by querying the DB per gene and ranking contexts by matched-gene count and average
  specificity.
- `PIASOmarkerDB(...)` — the underlying REST client class that the module-level functions
  wrap.

**Returned column schema:** `cell_type, condition, gene, species, specificity_score,
study_publication, tissue`.

```python
import piaso
studies = piaso.tl.queryPIASOmarkerDB(list_studies=True)      # list[str], 36 studies
df = piaso.tl.queryPIASOmarkerDB(cell_type="Microglia", limit=5)
# columns: cell_type, condition, gene, species, specificity_score, study_publication, tissue
df2 = piaso.tl.getMarkers(gene="AIF1", limit=3)               # alias of queryPIASOmarkerDB
res = piaso.tl.analyzeMarkers(["AIF1","P2RY12","CX3CR1","CSF1R"], n_top_genes=50)  # infer cell type
```

---

## Preprocessing utilities — `piaso.pp`

Besides `rotateSpatialCoordinates` (above), the preprocessing module provides small
data-wrangling helpers:
- `piaso.pp.table(values, rank=False, ascending=False, as_dataframe=False)` — R-style
  `table()`: value counts of a categorical, optionally sorted, returned as a dict or DataFrame.
- `piaso.pp.getCrossCategories(df, col1, col2, delimiter='@')` — build an ordered pandas
  Categorical from the cross-combination of two columns (e.g. `batch@celltype`), respecting
  existing category orders.

## Plotting — `piaso.pl`

PIASO ships plotting helpers (matplotlib/seaborn/scanpy-based — remember the `matplotlib<3.9`
pin from the install section):
- `plot_embeddings_split(adata, color, splitby, basis='X_umap', ...)` — faceted embedding
  scatter, one panel per `splitby` category, with shared axes and legend.
- `plot_features_violin(adata, feature_list, groupby=None, ...)` — stacked per-feature violins.
- `plotConfusionMatrix(data, groupby_query, groupby_reference, normalize='query', ...)` — an
  SVD-reordered confusion-matrix heatmap between two label columns.
- `createCustomCmapFromHex(hex_colors)` — build a colormap from hex colors; the
  `piaso.pl.color` submodule also exposes ready-made discrete/continuous palettes.
- `plotLigandReceptorInteraction` / `plotLigandReceptorLollipop` — consume the `runSCALAR`
  output (see the SCALAR section above).

> The `piaso.pp` / `piaso.pl` entries are documented from source; the plotting calls are not
> part of the executed test suite (they need a display and figure fixtures).

## See also

- `components/cosg.md` — the marker-gene method PIASO depends on (you still `import cosg`).
- `components/laris.md` — spatial ligand–receptor analysis (the LR counterpart to SCALAR).
- `gotchas.md` — ecosystem-wide traps (the `matplotlib<3.9` import bug, the annoy segfault,
  raw-vs-normalized layer requirements, PIASOmarkerDB needing internet).
