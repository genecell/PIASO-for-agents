# cytome — component reference (Python + R, self-sufficient)

**cytome** is the ecosystem's storage layer: one SQLite-backed **`.cytome` file** per dataset
holding expression matrices (chunked, compressed CSR), cell / gene / peak tables (SQL-queryable),
embeddings, kNN graphs, ATAC fragments (with an R*-tree range index), tissue images with
spatial-coordinate indexing, a JSON metadata store and a provenance log. It is *infrastructure*,
not a method: **every PIASO, COSG, LARIS and cytorete function accepts a `.cytome` path or open
Dataset wherever it accepts an AnnData**, reads it in chunks so peak memory is set by `batch_size`
rather than cell count (validated to 1.5 M cells), and writes results back onto the file. The
**R package** (`cytome`, r-universe) reads, writes and streams the same file into
SingleCellExperiment or Seurat with no Python. Format spec:
https://github.com/genecell/cytome/blob/main/docs/FORMAT_SPEC.md. Tutorials:
https://piaso.org/tutorials/cytome-basics/, `cytome-conversion`, `cytome-in-r`.

## Install

```bash
pip install cytome                 # Python — numpy/scipy/pandas/lz4 only (sqlite3 is stdlib); installed with piaso-tools anyway
pip install 'cytome[anndata]'      # + anndata for from_anndata / to_anndata
```
```r
install.packages("cytome", repos = "https://genecell.r-universe.dev")   # binaries (Windows/macOS); source build needs liblz4, libzstd, zlib
```

Executed against **cytome 0.3.1** (Python) and documented from **cytome (R) 0.1.0**. The R
package is MIT; the Python package BSD-3-Clause. No paper yet — cite the repositories.

## When to use a cytome (decision rule)

- **AnnData**: fits in RAM, one sample, you want scanpy interop.
- **`.cytome`**: more than ~100k cells; several samples in one object; ATAC fragments; a tissue
  image with region-of-interest queries; results that must persist on disk between sessions; or
  an **R collaborator** (Seurat / SCE ↔ Python with no bridge process).
- The analysis code is the same either way — only the object you pass changes.

## Naming contract (the thing to know before anything else)

| AnnData | cytome | note |
|---|---|---|
| `adata.X` (raw integers) / `layers["counts"]` | matrix **`RNA_counts`** | `{modality}_counts` holds **raw integers or does not exist** (≥ 0.3.0) |
| `adata.X` (normalized, non-integer) | matrix `RNA_data` (name recorded in `_anndata_X_layer`) | `cosg.cosg(path, layer="auto")` reads whichever the file records |
| `layers["infog"]` | `RNA_infog` | `piaso.tl.infog(ds, save_layer=True)` writes it |
| `obs` | `ds.cells` (SQL table; `.to_pandas()`, `.query("...")`, `.query_mask("...")`, `ds.cells["col"]`) | metrics, clusters, predictions land here |
| `var` | `ds.genes` (or `ds.peaks`) | |
| `obsm["X_umap"]`, `obsm["spatial"]`, `obsm["X_gdr"]` | embeddings **`RNA_umap`**, **`RNA_spatial`**, `RNA_gdr` | the `X_` / `obsm` tokens are dropped on conversion; `to_anndata` restores them; PIASO writes its own keys as given (`X_umap`, `X_gdr2`) and resolves either form |
| `obsp` graphs | `ds.graphs` | shared `graph_edges` table with R |
| `uns` | `ds.metadata` (JSON store) | e.g. cytorete's `regulon` entry |
| `uns["spatial"]` images | `ds.spatial_images` | `add_spatial_image`, `crop`, `as_uns` |

Files written before 0.2.6 use `RNA_obsm_X_umap`; both forms are readable
(`list(ds.embeddings.keys())` shows which).

## Python — building a cytome (executed)

```python
import cytome, piaso
# 1. from an AnnData already in memory
ds = cytome.from_anndata(adata, modality="RNA", output="e18.cytome")     # returns the OPEN dataset — keep using it
ds.n_cells, ds.n_genes, ds.list_matrices(), list(ds.embeddings.keys())
# (5865, 32285, ['RNA_counts', 'RNA_infog'], ['RNA_gdr', 'RNA_svd', 'RNA_umap', ...])

# 2. from an .h5ad on disk WITHOUT loading the matrix (the route for data that does not fit)
ds = cytome.from_h5ad("big.h5ad", output="big.cytome", modality="RNA", backed=True,
                      chunk_size=2048, storage_chunk_size=128)

# 3. from Cell Ranger / 10x output — no AnnData in between
ds = cytome.from_10x_h5("filtered_feature_bc_matrix.h5", "sample.cytome", sample_name="S1")   # RNA (+ ATAC for multiome)
# ds = cytome.from_cellranger("outs/", output="sample.cytome"); cytome.from_cellranger_arc("outs/", output=..., import_fragments=True)
# piaso.pp.importCellRanger(path, output=) — PIASO's Rust importer for RNA + ATAC fragments

# 4. several samples -> one file
merged = cytome.merge(["s1.cytome", "s2.cytome"], output="merged.cytome", batch_key="sample_id")
```

**One writer per file.** Use the Dataset the constructor returns; opening the same path a second
time while it is open leaves two writers and the next write fails (`database is locked`). Call
`ds.flush()` after writing columns/embeddings yourself and `ds.close()` when done (`cytome.open(path)`
reopens). `ds.copy()` / `ds.backup()` snapshot safely; `ds.checkpoint()` before `cp`/`rsync`.

## Python — reading and querying (executed)

```python
ds = cytome.open("e18.cytome")
ds.cells.columns                                   # ['cell_idx', 'barcode', 'sample_id', 'n_counts', 'leiden', 'pred', ...]
sub = ds.cells.query("leiden = '0'")               # SQL WHERE -> DataFrame (408 rows)
mask = ds.cells.query_mask("leiden IN ('0','1')")  # boolean mask for streaming / subsetting
block = ds.RNA.counts[:100, :50]                   # random access -> scipy CSR (reads only the needed chunks)
rows = ds.RNA.counts.rows(np.where(mask)[0])
for start, end, chunk in ds.RNA.counts.iter_rows():  # bounded-memory streaming, one CSR chunk at a time
    ...
umap = ds.embeddings["RNA_umap"]                   # numpy (n_cells, 2)
ds.set_categories("leiden", order=[...], colors={...})   # category order + palette persisted IN the file; every later plot honours it
ds.metadata.keys(); ds.provenance.show()           # JSON store; provenance log (parameters, versions, methods text)
```

## Python — analysing a cytome with the ecosystem (executed)

Same calls as on an AnnData; pass `modality="RNA"` on multimodal files.

```python
piaso.pp.calculateCellMetrics(ds, modality="RNA", prefix_vars={"mt": "mt-", "ribo": ["Rps", "Rpl"]})   # -> ds.cells columns
piaso.pp.scrublet(ds, expected_doublet_rate=0.06)
piaso.tl.infog(ds, modality="RNA", n_top_genes=3000, save_layer=True)     # save_layer=True writes RNA_infog (costs disk); False recomputes on the fly
piaso.tl.runSVD(ds, modality="RNA", layer="infog", n_components=50, key_added="X_svd")
piaso.tl.neighbors(ds, use_rep="X_svd"); piaso.tl.leiden(ds, resolution=1.0, key_added="leiden"); piaso.tl.umap(ds, use_rep="X_svd")
res = cosg.cosg(ds, groupby="leiden", modality="RNA", layer="infog", n_genes_user=25)   # RETURNS a dict: names, scores, groups_order
piaso.tl.runGDR(ds, groupby="leiden", layer="infog", modality="RNA", key_added="X_gdr")
piaso.tl.score(ds, gene_list=GABA, key_added="gaba", modality="RNA", layer="infog")    # -> ds.cells['gaba']
piaso.pl.embedding(ds, basis="X_umap", color="leiden"); piaso.pl.dotplot(ds, genes, groupby="leiden", modality="RNA", cytome_layer="infog")
piaso.pp.filter_cells(ds, mask=np.asarray(mask), inplace=False, output="subset.cytome")   # a cytome subset is a COPY; cell_mask= on plots selects without writing
adata_back = ds.to_anndata(modality="RNA")         # obsm keys and layers restored ('X_umap', 'X_gdr', layers['infog'])
```

Two things differ from the in-memory path: `cosg.cosg` **returns** its result (no `.uns`), with
columns in `groups_order` — reindex any per-cluster table by that order; and `layer="auto"` in
COSG means "normalize the stored counts as you stream", whereas a named layer is read as-is.
LARIS (`prepareLRInteraction(ds, ...)`, `runLARIS(..., data=ds)`) and cytorete
(`inferRegulon(ds, ...)`) take the same object. Tissue images and ROI queries:
`workflows/spatial_transcriptomics.md`.

## R — `cytome` (r-universe)

A **reader, writer and streaming interface**, deliberately not an analysis package: build,
merge, subset and filter on the Python side; analyse in Seurat / Bioconductor or hand back to
Python.

```r
library(cytome)
sce <- read_cytome("data.cytome")                      # SingleCellExperiment (default); RNA+ATAC -> ATAC altExp; embeddings -> reducedDims; graphs -> colPairs
so  <- read_cytome("data.cytome", as = "Seurat")       # Seurat (ATAC as a second assay; embeddings -> DimReduc; graphs slot)
x   <- read_cytome("data.cytome", as = "cytome")       # the open handle — look before loading
cytome_matrices(x); cytome_obs(x, "cells"); cytome_var(x, "RNA_counts"); cytome_embeddings(x); cytome_graphs(x)
M   <- read_cytome_matrix(x, "RNA_counts")             # features x cells dgCMatrix
cytome_close(x)

write_cytome(so,  "out.cytome")                        # generic: dispatches on Seurat or SCE; counts + embeddings + graphs by default
write_cytome(sce, "out.cytome", layers = TRUE)         # also carry normalized assays (off by default: recomputable, and floats compress badly)
write_cytome(so,  "out.cytome", graphs = FALSE)

# data too large to load
sce <- read_cytome("big.cytome", delayed = TRUE)       # DelayedArray-backed assay, chunk-aligned; scran/scuttle block-process it
totals <- cytome_stream(x, "RNA_counts", function(chunk, i0, i1, k) Matrix::rowSums(chunk))   # your own per-chunk reduction
```

What travels (measured on the reference file): counts and a second modality (`RNA` → `genes`,
`ATAC` → `peaks`; other modalities must be built in Python), embeddings, **graphs by default**,
**normalized layers opt-in**, all cell annotations; feature annotations beyond id + symbol do
**not** travel. ATAC feature ids must look like `chr1:100-200`. Cross-language conformance is
tested in CI in **both** directions on a Python-written reference file.

## The Seurat ↔ AnnData bridge (the conversion people actually want)

```r
write_cytome(seurat_obj, "shared.cytome")              # R
```
```python
adata = cytome.open("shared.cytome").to_anndata(modality="RNA")   # Python — or just pass the path to any piaso/cosg call
cytome.from_anndata(adata, output="shared.cytome")                 # back
```
```r
so <- read_cytome("shared.cytome", as = "Seurat")      # R
```

No `reticulate`, no HDF5 intermediate; either side can be the one that never opens the other
language. This is also how an R user reaches the Python-only methods (INFOG, GDR, PIASOscore,
SCALAR, LARIS, cytorete, COSG p-values).

## CLI

`cytome convert | info | merge | subset | downsample | export | validate | provenance | copy` —
`cytome info file.cytome` is the fastest way to see what a file holds.

## See also

- `workflows/streaming_large_data.md` — the end-to-end pipeline on a file, when to convert, R handoff.
- `workflows/spatial_transcriptomics.md` — `add_spatial_image`, `set_spatial_coords`, `cells_in_region`.
- `components/cosg.md` (dict return, `output_format`), `components/piaso.md` (`data=` contract).
- `gotchas.md` — counts invariant, embedding renames, one writer per file, `layer="auto"`.
