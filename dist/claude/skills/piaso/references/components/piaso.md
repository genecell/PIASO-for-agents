# PIASO — component reference

PIASO (`piaso-tools`) is the analysis layer of the ecosystem: a **self-contained Python + Rust
toolkit** for scRNA-seq, snRNA-seq and spatial transcriptomics that takes raw counts to figures
with **no scanpy dependency** — reading 10x data, QC, doublets, **INFOG** normalization, SVD and
**GDR** (marker-gene-guided dimensionality reduction), neighbors / Leiden / UMAP, gene-set
scoring (**PIASOscore**), marker- and reference-based cell-type prediction, **SCALAR**
single-cell ligand–receptor inference, a client for the remote **PIASOmarkerDB**, a plotting
suite, and `piaso.data` for datasets, genomes, motif databases and CellChatDB. Every function
accepts an **AnnData or a `.cytome`** (open Dataset or path) and streams on the latter. This file
is self-sufficient: it assumes nothing is already installed or imported. Tutorials:
https://piaso.org/tutorials/ · generated API reference: https://piaso.org/api/.

## Install

```bash
pip install piaso-tools          # also installs cosg and cytome (hard dependencies)
# optional extras
pip install 'piaso-tools[scanpy]'   # scanpy interop only — the core workflow needs none of it
pip install 'piaso-tools[harmony]'  # piaso.tl.runHarmony
```

Wheels ship for Linux, macOS (Intel + Apple Silicon) and Windows, Python 3.9–3.12; nothing to
build. **Do not pin matplotlib** — the `matplotlib<3.9` requirement applied to ≤ 1.1.0 only.
Blocks in this file were executed against **piaso-tools 1.2.3**; `predictCellTypeByGDR` needs
≥ 1.2.3 (earlier 1.2.x discarded its result on AnnData). Conda (`bioconda::piaso`) lags PyPI —
prefer pip.

## Import convention and the data= contract

```python
import piaso
piaso.tl   # tools          (== piaso.tools)      normalization, DR, clustering, scoring, annotation, SCALAR, PIASOmarkerDB
piaso.pp   # preprocessing  (== piaso.preprocessing)  reading, QC metrics, doublets, filtering, spatial coordinate helpers
piaso.pl   # plotting       (== piaso.plotting)   embeddings, dot plots, violins, scatter, sankey, stacked bars
piaso.data #                datasets, genome references, motif DBs, CellChatDB
piaso.settings.set_figure_params(style="cell")   # one house style; 'nature', 'science', ... also exist
```

Short and long forms are the **same runtime objects**. The first argument of every `pp` / `tl` /
`pl` function is `data` (the old `adata=` still works) and may be an `AnnData`, an open
`cytome.CytomeDataset`, or a path to a `.cytome` file; on a cytome the function streams in
`batch_size` chunks and writes results back onto the file (`cells` table, embeddings,
`metadata`). Pass `modality="RNA"` on multi-modal cytomes. See `components/cytome.md`.

## Dependency: COSG

`piaso-tools` hard-depends on `cosg` and calls `cosg.cosg(...)` inside GDR, `specificity_matrix`,
`cospecificity_trans` and `stitchSpace` — but COSG is **not re-exported** under `piaso.*`. To call
it yourself, `import cosg`. See `components/cosg.md`.

## Citation

> Wu, S.J., Dai, M. et al. Pyramidal neurons proportionately alter cortical interneuron
> subtypes. Nature (2026). DOI: 10.1038/s41586-025-09996-8

---

## Shared setup for the examples (executed)

All blocks below continue from this state — the `e18_v3_nuclei` fixture (10x h5, ~6k mouse brain
nuclei, see `data.md`) read with PIASO's own reader. `.X` holds **raw UMI counts**, which is what
`infog` expects.

```python
import numpy as np, pandas as pd
import piaso, cosg
adata = piaso.pp.read_10x_h5("e18_v3_nuclei.h5")     # or: piaso.data.load_dataset("e18_v3_nuclei")
# (5973, 32285); var columns gene_ids / gene_symbols / feature_types / genome; names made unique
```

---

## Reading data and QC — `piaso.pp`

- `piaso.pp.read_10x_h5(path, modality="rna", var_names="gene_symbols", make_unique=True)` — Cell
  Ranger HDF5 → AnnData with raw counts in `.X`. `piaso.pp.read_10x(path)` accepts an `.h5` or an
  MTX directory; `read_10x_mtx` for the directory form. `piaso.pp.importCellRanger(path, output=)`
  builds a **cytome** (RNA + ATAC fragments, Rust import) instead of an AnnData.
- `piaso.pp.calculateCellMetrics(data, prefix_vars={...}, feature_set_vars=None)` — writes
  `obs["n_counts", "n_genes"]` and, per prefix group, `n_counts_<k>` / `pct_counts_<k>`.
  **Prefixes are species-cased**: mouse `"mt-"`, `["Rps", "Rpl"]`; human `"MT-"`, `["RPS", "RPL"]`.
  Matching nothing yields zeros silently.
- `piaso.pp.scrublet(data, expected_doublet_rate=0.06, library_key=None, random_state=0)` — PIASO's
  own Scrublet: writes `obs["scrublet_score", "is_doublet"]`. With several libraries pass
  `library_key=` so each is scored against its own background. The automatic threshold is
  conservative; a fixed cut (`scrublet_score > 0.3`) is common practice.
- `piaso.pp.filter_cells(data, min_counts=, max_counts=, min_features=, max_features=, mask=)` —
  in place on AnnData; on a cytome `inplace=False, output=` writes a new file (a cytome subset is
  a copy). `filter_features` for the gene side; `calculateGroupMetrics(data, groupby=)` for a
  per-group QC table (feed to `piaso.pl.plotGroupMetrics`).

```python
piaso.pp.calculateCellMetrics(adata, prefix_vars={"mt": "mt-", "ribo": ["Rps", "Rpl"]})
piaso.pp.scrublet(adata, expected_doublet_rate=0.06, random_state=0)
adata = adata[~adata.obs["is_doublet"].astype(bool)].copy()
piaso.pp.filter_cells(adata, min_counts=500, min_features=250)
# obs: n_counts, n_genes, n_counts_mt, pct_counts_mt, n_counts_ribo, pct_counts_ribo, scrublet_score, is_doublet
```

Filter **both tails** (`n_counts` above the 99th percentile is often doublets) and read the
distribution before choosing a mitochondrial cut: nuclei samples sit near 0 %, whole cells near
10 %+. After removing cells or clusters, **re-run INFOG → SVD → neighbors → Leiden** — the gene
selection and SVD were fitted with the removed cells present.

---

## INFOG normalization — `piaso.tl.infog`

**What it computes:** information-content normalization of **raw integer UMI counts** plus
selection of the top-`n_top_genes` informative genes, in one step. **Reads** `adata.X` or
`layer=` (refuses non-integer input and names the layer to pass; `allow_non_integer=True`
overrides). **Writes** `adata.layers["infog"]` (or `key_added`), `adata.var["infog_var"]`,
`adata.var["highly_variable"]`. On a cytome: `save_layer=True` persists `RNA_infog` on the file,
otherwise it is recomputed on the fly by later steps.

```python
piaso.tl.infog(adata, n_top_genes=3000)                 # counts in .X (10x h5)
# piaso.tl.infog(adata, layer="raw", n_top_genes=3000)  # when the registry says counts_layer='raw'
```

Datasets from the registry may keep counts in a layer — check
`piaso.data.dataset_info(name)["counts_layer"]` (`sea_ad_mtg_20k` → `UMIs`,
`adult_cortex_multiome_rna` → `raw`). `piaso.tl.infog_svd(...)` does INFOG + HVG + SVD in one call;
the old `runSVDLazy` is a deprecated alias of it.

---

## Embedding and clustering — `runSVD`, `neighbors`, `leiden`, `umap`

Standard steps, PIASO-native (igraph Leiden; no scanpy, no leidenalg needed):

```python
piaso.tl.runSVD(adata, layer="infog", n_components=50, key_added="X_svd")   # ALWAYS pass layer= — default is .X (raw counts)
piaso.tl.neighbors(adata, use_rep="X_svd", n_neighbors=15)                  # writes obsp['distances','connectivities'], uns['neighbors']
piaso.tl.leiden(adata, resolution=1.0, key_added="leiden")                  # writes obs['leiden'] (categorical) — lowercase by default
piaso.tl.umap(adata, use_rep="X_svd")                                       # writes obsm['X_umap']
# 22 clusters on the fixture at resolution 1.0
```

`neighbors(key_added=)` / `leiden(neighbors_key=)` / `umap(neighbors_key=, key_added=)` let you
keep a second graph (e.g. on `X_gdr`) beside the first. `piaso.tl.runHarmony(adata, batch_key=,
use_rep="X_svd", key_added=)` corrects an embedding when a real batch effect exists (measure it
first — see the multi-sample tutorial; GDR with `batch_key` is usually the better route).

---

## Marker-gene-guided DR (GDR) — `piaso.tl.runGDR` / `runGDRParallel` / `projectGDR`

**What it computes:** GDR is **marker Gene-guided Dimensionality Reduction**. It (1) takes groups
(`groupby`, or clusters within each `batch_key` when `groupby=None`), (2) runs **COSG** for the
top-`n_gene` markers per group, (3) scores every cell against each marker set with PIASOscore,
(4) L2-normalizes. **That cell × marker-group score matrix is the embedding** — its width is the
number of groups, not `n_svd_dims`. With `batch_key`, markers are selected **within each batch**
and every cell is scored against all batches' sets, which integrates libraries by identity
rather than by correction.

**Reads** `layers[layer]` (normalized; default `infog`), `obs[groupby]` / `obs[batch_key]`; with
`batch_key`, `infog_layer=<counts layer>` when `.X` is not raw counts (GDR re-runs INFOG per
batch — pointing it at a scaled `.X` produces NaN). **Writes** `obsm[key_added]` (default `X_gdr`),
`uns["gdr"]`, and with `save_reference=True` (default) `uns["gdr_reference"]` — the marker sets +
column norms that define the space. Defaults: `n_gene=20`, `mu=10.0`.

```python
piaso.tl.runGDR(adata, groupby="leiden", layer="infog", key_added="X_gdr")    # (n_obs, 22) on the fixture
piaso.tl.neighbors(adata, use_rep="X_gdr", n_neighbors=15, key_added="gdr")
piaso.tl.umap(adata, use_rep="X_gdr", key_added="X_umap_gdr", neighbors_key="gdr")
# multi-sample: piaso.tl.runGDR(adata, batch_key="sample", groupby=None, layer="infog", n_gene=30, key_added="X_gdr")
```

`runGDRParallel` is the multi-process variant with the same contract (`max_workers`,
`calculate_score_multiBatch`). **`projectGDR`** freezes a reference's space and projects new data
into it — no re-fitting, reference coordinates stay fixed:

```python
query = adata[np.random.default_rng(0).random(adata.n_obs) < 0.3].copy()      # any AnnData / .cytome with overlapping var_names
piaso.tl.projectGDR(query, reference=adata, key_added="X_gdr_proj")          # needs uns['gdr_reference'] on the reference
# query.obsm['X_gdr_proj'] has the REFERENCE's width; mode='reference' (comparable coords) | 'self' (new-batch scaling)
```

Then `predictCellTypeByMarker(query, ..., use_rep="X_gdr_proj")` for label transfer. On a `.cytome`
query the projection streams. Median cosine to a jointly computed embedding is ~0.97 on the
cortex tutorial; if you need exact joint geometry, embed jointly.

---

## Gene-set scoring (PIASOscore) — `piaso.tl.score`

**What it computes:** per-cell enrichment of a gene set against **expression-matched control
sets** (each query gene's k nearest neighbours in mean/variance space; `n_ctrl_set=100`),
`score = score_query − score_ctrl_average`, with Monte-Carlo and pooled empirical p-values + BH
FDR when `compute_pvalues=True`. All sets and their controls are evaluated as one sparse matmul
per chunk in a **Rust kernel**, so a 300-set pathway database costs about the same as one set.
Removes the depth correlation a naive mean carries (−0.24 → 0.01 on SEA-AD).

**Reads** `layers[layer]` (default `infog` — errors if absent). **Single set** → writes
`obs[key_added]` + `uns[key_added]` (11 columns with p-values) and returns `None`. **Dict /
DataFrame / list-of-lists** → **returns a tuple** `(score_matrix, names, pvals)`; `compute_pvalues`
is off by default there.

```python
GABA = [g for g in ["Gad1", "Gad2", "Slc32a1", "Dlx1", "Dlx5"] if g in adata.var_names]
piaso.tl.score(adata, gene_list=GABA, key_added="gaba", compute_pvalues=True)
# obs['gaba']; uns['gaba'].columns == ['score','score_query','score_ctrl_average','pval_mc','nlog10_pval_mc',
#   'pval_mc_FDR','nlog10_pval_mc_FDR','pval','nlog10_pval','pval_FDR','nlog10_pval_FDR']
sig = adata.uns["gaba"]["pval_mc"] < 0.01            # a per-cell gate

sets = {"gaba": GABA, "glut": [g for g in ["Slc17a7", "Neurod6", "Satb2"] if g in adata.var_names]}
score_matrix, names, pvals = piaso.tl.score(adata, gene_list=sets, layer="infog")   # (n_obs, 2), list, None
```

The score matrix is cells × sets — the same shape as expression, so `anndata.AnnData(X=score_matrix,
obs=adata.obs, var=pd.DataFrame(index=names))` + `cosg.cosg(..., groupby=)` tells you **which
pathway is a marker of which cell type** (KEGG via `gseapy.parser.get_library`, drug targets via
`piaso.data.load_chembl_targets()`; case-fold symbols for mouse). On a cytome the same call streams.

---

## Cell-type prediction — `predictCellTypeByMarker` / `predictCellTypeByGDR`

### `predictCellTypeByMarker` — from marker sets

Scores every cell against each cell type's marker set (PIASOscore on `score_layer`), assigns the
argmax, then smooths over the kNN graph of `use_rep` (default `X_gdr` — **pass `use_rep="X_svd"` if
that is all you have**). Marker sets come from COSG on your own clusters or a reference, or from
PIASOmarkerDB. **Writes** `obs[key_added]` (default `CellTypes_predicted`) plus `_raw`,
`_smoothed`, `_score`, `_smoothed_confidence`, `_confidence_smoothed` columns and
`obsm[key_added + "_score"]`. Labels are the **keys** of `marker_gene_set`.

```python
# from a curated study — as_dict=True returns a TUPLE (table, dict): unpack both
markers_df, marker_sets = piaso.tl.getMarkers(study="AllenWholeMouseBrain_isocortex", as_dict=True)   # 1300 rows, 26 types
piaso.tl.predictCellTypeByMarker(adata, marker_gene_set=marker_sets, score_layer="infog",
                                 use_rep="X_gdr", key_added="CellTypes")
adata.obs["CellTypes"].value_counts().head()        # '001 CLA-EPd-CTX Car3 Glut', '005 L5 IT CTX Glut', '319 Astro-TE NN', ...
```

Labels keep the source taxonomy's names, which is a feature: they are traceable to the study.
Check the result against its own markers with `piaso.pl.dotplot(adata, <top genes>,
groupby="CellTypes", standard_scale="var")`. Held-out accuracy on the cortex tutorial: 95.4 %
with COSG `mu=10` markers on a GDR embedding.

### `predictCellTypeByGDR` — reference-based label transfer (≥ 1.2.3)

Concatenates reference and query, runs GDR jointly, Harmony-integrates `X_gdr`, trains an SVM
on the reference labels and predicts the query. **Reads** `adata.layers[layer]` /
`adata_ref.layers[layer_reference]` (default `"log1p"` — pass `"infog"` if that is what you have),
`adata_ref.obs[reference_groupby]`, `adata.obs[query_groupby]`. **Writes**
`adata.obs[key_added or "CellTypes_gdr"]`. Two contracts found by running it: query and reference need **disjoint cells
with distinct `obs_names`**, and the `reference_groupby` column must **also exist in the query** (a placeholder
value is fine) — otherwise `KeyError: '<reference_groupby>'`.

```python
query.obs["CellTypes"] = "unknown"                                    # placeholder: the reference_groupby column must exist on the query
piaso.tl.predictCellTypeByGDR(query, adata_ref, layer="infog", layer_reference="infog",
                              reference_groupby="CellTypes", query_groupby="leiden", key_added="CellTypes_gdr")
```

Prefer `runGDR(save_reference=True)` + `projectGDR` + `predictCellTypeByMarker` when the
reference's coordinates must stay fixed across many queries.

---

## Local sub-clustering — `piaso.tl.leiden_local`

Re-clusters **selected groups** after recomputing the embedding on just those cells (unlike
scanpy's `restrict_to`, which reuses the global graph). **Writes** `obs[key_added]` (default
`Leiden_local`) with `{group}-{local}` labels.

```python
piaso.tl.leiden_local(adata, groupby="leiden", groups=["0"], resolution=0.2,
                      key_added="leiden_local", dr_method="X_svd")
```

`dr_method` accepts `X_svd`, `X_svd_full`, `X_pca`, GDR / Harmony variants (the latter need
`batch_key`).

---

## Marker-guided batch integration — `piaso.tl.stitchSpace`

**Not spatial** despite the name ("Space" = embedding space). Builds a BBKNN graph across
batches on `use_rep`, prunes edges between clusters whose COSG markers do not overlap, and applies
one centroid correction to the embedding. **Writes** `obsm[key_added]` (default `X_stitch`).
Requires `bbknn`; pin `annoy==1.16.3` (BBKNN segfaults with annoy ≥ 1.17). For most multi-sample
data `runGDR(batch_key=)` is the simpler route.

---

## SCALAR — single-cell ligand–receptor (`piaso.tl.runSCALAR`)

**SCALAR is a PIASO function, not a package.** It infers cell-type-resolved ligand–receptor
interactions for **dissociated** data: `interaction_score = specificity[ligand, sender] ×
specificity[receptor, receiver]`, with a permutation null built from expression-matched control
genes (`n_nearest_neighbors=30` matched genes, `n_permutations=1000` by default) and BH FDR **per
sender–receiver pair**. Two inputs, both now provided by the ecosystem:

- **`specificity_matrix`** (genes × cell types): `piaso.tl.specificity_matrix(data, groupby=,
  cosg_layer="counts")` runs COSG at full `n_genes_user` and pivots it into a dense frame (needs
  the raw-counts layer named by `cosg_layer`; on a cytome it reuses the cached COSG run).
  Equivalent by hand: `cosg.cosg(adata, groupby=, n_genes_user=adata.n_vars, mu=1,
  remove_lowly_expressed=False)` then pivot `uns['cosg']`. Scores are comparable **within** a
  sender–receiver pair, not across pairs; whether to rescale columns is an open question on the
  maintainers' side — this hub documents the raw COSG matrix.
- **`lr_pairs`**: `piaso.data.load_lr_database("mouse" | "human", annotation=None)` fetches
  **CellChatDB** (mouse 3105 / human 2951 pairs; columns include `ligand`, `receptor`,
  `pathway_name`, `annotation` ∈ {Secreted Signaling, ECM-Receptor, Cell-Cell Contact, Non-protein
  Signaling}) — the same tables LARIS bundles. Any DataFrame with ligand/receptor columns works.

**Reads** `layers[layer]` for the background statistics. **Writes nothing to `adata`** — returns a
DataFrame: `ligand, receptor, sender, receiver, interaction_score, p_value, p_value_fdr,
nlog10_p_value_fdr` (+ `annotation_col` if given).

```python
spec = piaso.tl.specificity_matrix(adata, groupby="leiden", cosg_layer="counts")   # (n_genes, n_groups); needs layers['counts']
lr = piaso.data.load_lr_database("mouse")                                          # (3105, 28) CellChatDB
res = piaso.tl.runSCALAR(adata, specificity_matrix=spec, lr_pairs=lr, layer="infog",
                         annotation_col="annotation", n_permutations=1000, random_seed=42)
sig = res[res["p_value_fdr"] < 0.05]
sig.groupby(["sender", "receiver"]).size().sort_values(ascending=False).head()
```

(On the ~6k-cell fixture with cluster labels nothing reaches FDR < 0.05; on the 20k-nucleus
SEA-AD cortex the tutorial reports 1.67 M tested interactions and VLMC → astrocyte collagen
signalling on top.) Interaction **counts track cell-type abundance and marker sharpness** —
compare pairs of comparable size. Plot one pair with the two specificities that produced it:

```python
sig = sig.copy()
sig["CellTypeXCellType"] = sig["sender"] + "@" + sig["receiver"]
sig["ligandXreceptor"] = sig["ligand"] + "-->" + sig["receptor"]
sig["ligand_specificity"] = [spec.at[r.ligand, r.sender] for r in sig.itertuples()]
sig["receptor_specificity"] = [spec.at[r.receptor, r.receiver] for r in sig.itertuples()]
piaso.pl.plotLigandReceptorInteraction(interactions_df=sig, specificity_df=spec,
                                       cell_type_pairs=["VLMC@Astrocyte"], ligand_receptor_sep="-->",
                                       top_n=30, y_max=float(sig["interaction_score"].max() * 1.1))
# piaso.pl.plotLigandReceptorLollipop(sig, cell_type_pairs=[...], col_cell_type_pair="CellTypeXCellType", ...)
```

**Decision rule — SCALAR vs LARIS:** spatial coordinates → LARIS (`components/laris.md`); none →
SCALAR. Same database either way; the natural pairing is SCALAR on the dissociated reference,
LARIS on the section. On a targeted panel count complete pairs first (both genes measured).

---

## PIASOmarkerDB — remote marker database client

**A REST API client, not bundled data**: HTTP calls to `https://piaso.org/piasomarkerdb/api/v1/`
(needs internet; caches under `~/.piaso/markers`). Python-only; from R or another agent use the
`piaso-mcp` server's `query_marker_db` or the API directly. 36 studies, human + mouse, brain,
blood, bone marrow, breast, thymus, whole-body atlases.

- `piaso.tl.getMarkers(gene=, cell_type=, study=, species=, tissue=, condition=, min_score=,
  max_score=, limit=, as_dict=False, list_studies=False, list_cell_types=False, list_genes=False)`
  → DataFrame with `cell_type, condition, gene, species, specificity_score (a COSG score),
  study_publication, tissue`; `as_dict=True` → **tuple** `(df, {cell_type: [genes]})`; the `list_*`
  flags → `list[str]`. `queryPIASOmarkerDB` is the same function.
- `piaso.tl.analyzeMarkers(genes, n_top_genes=50, species=, tissue=, studies=, exclude_studies=,
  exclude_cell_types=)` — the reverse query: a gene list → ranked DataFrame of matching cell types
  (`cell_type, study_publication, species, tissue, condition, matched_gene_count, matched_genes,
  avg_specificity`). A `{cluster: [genes]}` dict or COSG DataFrame → **tuple** `(results_dict,
  top_hits)` with `top_hits[cluster]` a cell-type **string** (`"Unassigned"` if none).

```python
studies = piaso.tl.getMarkers(list_studies=True)                                     # 36 study keys
df = piaso.tl.getMarkers(gene="Sst")                                                 # where is Sst a marker, across studies
markers_df, marker_sets = piaso.tl.getMarkers(study="AllenWholeMouseBrain_isocortex", as_dict=True)
names = pd.DataFrame(adata.uns["cosg"]["names"])                                     # after cosg.cosg(adata, groupby="leiden", ...)
results, top_hits = piaso.tl.analyzeMarkers({c: list(names[c]) for c in names.columns}, species="Mouse")
adata.obs["celltype_db"] = adata.obs["leiden"].map(top_hits)                        # e.g. '0' -> '038 DG-PIR Ex IMN'
```

Cell-type names must match the study's vocabulary exactly (`list_cell_types=True` first).

---

## Reference data — `piaso.data`

```python
piaso.data.list_datasets(); piaso.data.dataset_info("sea_ad_mtg_20k")["counts_layer"]   # 'UMIs'
adata = piaso.data.load_dataset("mouse_brain_10k_gemx")                                 # download, md5-verify, open
ds = piaso.data.load_dataset("sea_ad_mtg_20k_cytome", return_type="cytome")             # open a .cytome atlas
path = piaso.data.fetch_dataset("e18_v3_nuclei")                                        # download only
piaso.data.fetch_genome("hg38"); piaso.data.fetch_2bit("hg38"); piaso.data.fetch_jaspar()  # cytorete / motif inputs
lr = piaso.data.load_lr_database("human", annotation="Secreted Signaling")               # CellChatDB slice
```

Cache root `~/.piaso/data/` (override `data_dir=` / `piaso.settings.data_dir` / `PIASO_DATA_DIR`);
registry refreshed every 24 h (`piaso.data.refresh_registry()`). Also: `load_chembl_targets`,
`load_cisbp`, `load_cistarget_motifs`, `fetch_screen`, `load_jaspar_meme`, `build_tf_motif_map`,
`extract_sequences`, and `piaso.pp.scan_motifs(pwms, sequences)` (Rust PWM scanner) — the motif
engine cytorete builds on. See `data.md`.

---

## Spatial helpers

PIASO has no spatial *statistics* (that is LARIS); it has coordinate and image plumbing that
works on AnnData or cytome:
- `piaso.pp.rotateSpatialCoordinates(data, angle_degrees, spatial_key="X_spatial", backup_spatial_key=)`
  and `piaso.pp.alignSpatialCoordinates(data, groupby=, spatial_key="spatial", with_std=False)`
  (centre each section on its centroid for split plots).
- `piaso.pl.plotEmbedding(ds, color=, basis="spatial", image=True, img_key=, cell_mask=)` draws a
  registered tissue image under the cells when the **cytome** holds one (`ds.add_spatial_image`);
  `ds.cells_in_region(x=, y=)` selects a rectangle. See `workflows/spatial_transcriptomics.md`.

---

## Plotting — `piaso.pl`

All take AnnData or cytome, `save=` a path, `show=False` for headless use; set the house style
once with `piaso.settings.set_figure_params(style="cell")`.

| Call | Use |
|---|---|
| `embedding(data, basis="X_umap", color=, legend_loc="right"\|"on_data"\|"both", palette=, cmap=, vmin_pct=, vmax_pct=, ncol=)` (= `plotEmbedding`) | any embedding, categorical or continuous, one or several colours; `image=True` on spatial cytomes |
| `plot_embeddings_split(data, color=, splitby=, basis=)` | one panel per category of `splitby` |
| `dotplot(data, features, groupby=, standard_scale="var", layer=)` | marker check — the diagonal is the point |
| `violin(data, features, groupby=)` / `plot_features_violin(data, features, groupby=)` | per-group QC and expression |
| `scatter(data, x=, y=, color=, logx=, logy=, marginals=True)` | QC scatter with marginals |
| `stackedBarplot(data, groupby=, splitby=)` | composition per group |
| `sankey(data, left=, right=, right_order=)` / `plotConfusionMatrix(data, groupby_query=, groupby_reference=)` | compare two labellings |
| `plotGroupMetrics(table, data=, groupby=)` | render a `calculateGroupMetrics` table |
| `plotLigandReceptorInteraction`, `plotLigandReceptorLollipop` | SCALAR results |
| `piaso.pl.color.d_color3 / d_color4 / ...`, `createCustomCmapFromHex([...])` | palettes |

Continuous colour picks PIASO's sequential map for a **gene** and `Spectral_r` for a numeric
**cell column** (so metadata never reads as expression); an explicit `cmap=` wins. `palette=` may
be a colormap name for an ordered categorical (ages, stages).

---

## Not documented here

`piaso.tl.runATACLazy`, `run_TFIDF`, `picco`, `quantify_peaks` and the rest of the scATAC surface
exist in the wheel but are **not released** (not on piaso.org's API pages; "coming" on the module
grid). `piaso.tl.inferRegulon` / `inferGRN` / `inferTFActivity` / `regulonActivity` /
`regulonSpecificity` are **forwarders to cytorete** — `import cytorete` (`components/cytorete.md`).

## See also

- `components/cosg.md` — the marker method every GDR / SCALAR / cytorete step calls.
- `components/cytome.md` — the file format every function above streams from.
- `workflows/end_to_end_scrnaseq.md` — the executed pipeline these blocks come from.
- `gotchas.md` — layer contracts, deprecated names, the `as_dict` tuple, species-cased prefixes.
