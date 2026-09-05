# Workflow — spatial transcriptomics (Xenium / Visium HD / MERFISH / Stereo-seq / Slide-tags)

From a spatial section to annotated clusters on the tissue image, then the downstream questions
spatial data exists for — **ligand–receptor where cells touch (LARIS)** and **regulons (cytorete)**
— with the coverage check that decides which of them a targeted panel can support. The pipeline
is the scRNA-seq one (`end_to_end_scrnaseq.md`) plus coordinates, and the natural container is a
**cytome** (coordinates indexed, image stored, everything streamed). Blocks marked *(executed)* ran
on the LARIS Slide-tags tonsil object (`adata_tonsil.h5ad`, 5,695 cells, `obsm['X_spatial']`,
`obs['cell_type']`; see `data.md`); blocks marked *(from tutorial)* are lifted verbatim from the
executed piaso.org tutorials *Xenium with tissue-image overlay*, *Xenium into a cytome directly*,
*Downstream in situ*, *Stereo-seq whole embryo* and *MERFISH sections*, whose data is downloaded
from the original providers and was not re-run here.

## Install

```bash
pip install piaso-tools laris cytorete    # + tifffile for OME-TIFF images
```

## Decision rules first

- **Coordinates → LARIS** for ligand–receptor (`ligand_receptor.md`); SCALAR is for dissociated data.
- **Panel coverage before analysis**: an LR pair needs *both* genes measured; a pathway score
  enough of the pathway; a regulon the TF *and* its targets. 313-gene Xenium breast panel:
  11 / 2951 LR pairs complete (0.4 %) — not an option; 5,006-gene Xenium Prime: 1,637 (53 %);
  whole-transcriptome in situ: 2,831 (96 %).
- **AnnData vs cytome**: a section of > 100k cells, a morphology image, ROI queries, or several
  sections in one object → cytome. Under that, an AnnData with `obsm['spatial']` is fine.
- **Coordinate key**: PIASO plots read `basis="spatial"` (cytome embedding `RNA_spatial` /
  `spatial`); LARIS reads `use_rep_spatial="X_spatial"`. Copy one to the other when mixing:
  `adata.obsm["X_spatial"] = np.asarray(adata.obsm["spatial"], dtype=np.float64)`.

## Step 0 — Coverage check *(executed on the tonsil; tutorial numbers above)*

```python
import numpy as np, pandas as pd
import laris
lr = laris.datasets.lrDatabase(species="human")                    # CellChatDB, 2951 pairs (== piaso.data.load_lr_database("human"))
pairs = lr[["ligand", "receptor"]].dropna().drop_duplicates()
genes = set(adata.var_names)
complete = pairs.apply(lambda r: r.ligand in genes and all(x.strip() in genes for x in str(r.receptor).split(",")), axis=1)
complete.sum(), len(pairs)                                          # tonsil (whole transcriptome): 1985 pairs have both genes measured
```

## Step 1 — Cells, coordinates, clustering *(from tutorial: Xenium breast, 164k cells)*

```python
import piaso, cosg
adata = piaso.pp.read_10x_h5("Xenium_..._cell_feature_matrix.h5")          # Xenium's matrix is 10x format
adata.var_names_make_unique()
cells = pd.read_csv("Xenium_..._cells.csv.gz")
cells["cell_id"] = cells["cell_id"].astype(str)                            # h5 obs_names are strings
cells = cells.set_index("cell_id").loc[adata.obs_names]                    # ALIGN on the barcode, never zip positionally
adata.obsm["spatial"] = cells[["x_centroid", "y_centroid"]].to_numpy()     # microns
adata = adata[np.asarray(adata.X.sum(1)).ravel() >= 10].copy()             # light QC
piaso.tl.infog(adata, n_top_genes=2000)
piaso.tl.runSVD(adata, layer="infog", n_components=30, key_added="X_svd")
piaso.tl.neighbors(adata, use_rep="X_svd", n_neighbors=15)
piaso.tl.leiden(adata, resolution=1.0)                                     # obs['leiden']
```
Or straight into a cytome without an AnnData *(from tutorial)*: `ds = cytome.from_10x_h5(
"cell_feature_matrix.h5", output="xenium.cytome")` (blank / negative-control codewords are
skipped with a warning — keep them for QC from the h5, not in the RNA modality), then
`ds.add_embedding("spatial", xy); ds.set_spatial_coords(xy); ds.flush()` with `xy` reindexed on
`ds.cells["barcode"]`, and the same `infog` / `runSVD` / `neighbors` / `leiden` calls on `ds`.

## Step 2 — One file: matrix + coordinates + index + image *(from tutorial)*

```python
import tifffile, cytome
with tifffile.TiffFile("..._morphology_mip.ome.tif") as tf:               # pyramidal OME-TIFF; one mid level (<= 8000 px) is plenty
    level = next(i for i, l in enumerate(tf.series[0].levels) if max(l.shape[:2]) <= 8000)
    img = tf.series[0].levels[level].asarray(); full_w = tf.series[0].levels[0].shape[1]
scalef = (img.shape[1] / full_w) / 0.2125                                  # micron -> stored pixel (Xenium full-res is 0.2125 um/px)
img8 = np.clip(img.astype(np.float32) / np.percentile(img, 99) * 255, 0, 255).astype(np.uint8)

ds = cytome.from_anndata(adata, output="xenium.cytome")                    # stores obsm['spatial'] as embedding + builds the R*-tree index
ds.add_spatial_image("xenium_rep1", "morphology", img8,
                     scalefactors={"tissue_morphology_scalef": scalef, "spot_diameter_fullres": 10.0})
# a PNG/JPEG/TIFF path works too: ds.add_spatial_image("xenium", "morphology", "morphology_mip.ome.tif")
```
The one number that must be right is the **scale factor**. The image is stored losslessly inside
the same SQLite file and travels with the data.

## Step 3 — Clusters and genes on the tissue *(from tutorial)*

```python
piaso.pl.plotEmbedding(ds, color="leiden", basis="spatial", image=True, img_key="morphology",
                       point_size=0.3, alpha=0.7, legend_loc="right")           # orientation/units handled: image drawn in coordinate space
piaso.pl.plot_embeddings_split(ds, color="leiden", splitby="leiden", basis="spatial",
                               image=True, img_key="morphology", ncol=5)        # one cluster per panel: which clusters have a PLACE
piaso.pl.plotEmbedding(ds, color="KRT14", basis="spatial", image=True, img_key="morphology", cmap="Spectral_r")   # a gene, same call
```

## Step 4 — Name the clusters, and a GDR embedding from their markers *(from tutorial)*

```python
markers = cosg.cosg(ds, groupby="leiden", modality="RNA", layer="infog", n_genes_user=5, mu=10, output_format="dict")   # dict return on a cytome
piaso.tl.runGDR(ds, batch_key=None, groupby="leiden", n_gene=20, mu=10, layer="infog", score_layer="infog", key_added="X_gdr")
piaso.tl.neighbors(ds, use_rep="X_gdr", n_neighbors=15, key_added="gdr")
piaso.tl.umap(ds, use_rep="X_gdr", key_added="X_umap_gdr", neighbors_key="gdr")
```
Name them with `piaso.tl.analyzeMarkers(...)` or a matched PIASOmarkerDB study
(`markerdb_annotation.md`); the tutorial reads KRT14/KRT5/MYLK → myoepithelium, CD14/MRC1/CD163 →
macrophages, MZB1/TNFRSF17 → plasma cells off the top markers.

## Step 5 — Regions of interest: cells and pixels from the same rectangle *(from tutorial)*

```python
xy = ds.embeddings["RNA_spatial"]
cx, cy = xy[ds.cells.to_pandas()["leiden"].astype(str).to_numpy() == "9"].mean(axis=0)
cells_in = ds.cells_in_region(x=(cx - 250, cx + 250), y=(cy - 250, cy + 250))          # indexed R*-tree lookup, a 500 um window
piaso.pl.plotEmbedding(ds, color="leiden", basis="spatial", image=True, img_key="morphology",
                       cell_mask=cells_in, point_size=6, alpha=0.85)                 # the query straight in; spatial_images.crop(...) cuts the matching pixels
```

## Step 6 — Several sections; orienting a section *(from tutorial)*

```python
piaso.pl.plot_embeddings_split(ds, color="Cell_class", splitby="section", basis="spatial", ncol=3)   # MERFISH: sections in one cytome
piaso.pp.rotateSpatialCoordinates(ds, angle_degrees=180, spatial_key="spatial", backup_spatial_key="spatial_delivered")   # Stereo-seq: orient
piaso.pp.alignSpatialCoordinates(ds, groupby="section", spatial_key="spatial", key_added="spatial_aligned")            # centre each section for split plots
```
Bin-resolution data (Stereo-seq, Visium HD) flows through the same steps; `plotEmbedding(...,
cell_mask=)` zooms by coordinates.

## Step 7 — Downstream: ligand–receptor where cells touch (LARIS) *(executed on the tonsil)*

```python
import laris
lr_df = laris.datasets.lrDatabase(species="human")
lr_data = laris.tl.prepareLRInteraction(adata, lr_df, use_rep_spatial="X_spatial")            # (cells x LR pairs) diffused over the spatial kNN graph
bg = laris.tl.prepareLRBackground(adata, lr_df, use_rep_spatial="X_spatial")                 # exact matched-gene null; slow, reusable — build once
laris_lr, res = laris.tl.runLARIS(lr_data, adata, use_rep="X_spatial", use_rep_spatial="X_spatial",
                                  groupby="cell_type", background=bg)
laris.pl.plotCCCSpatial(lr_data, basis="X_spatial", interaction=res.iloc[0]["interaction_name"], color_by="score")   # the score, on the section
```
Group by **clusters**, not hundreds of fine types (LARIS tests every sender–receiver pair). Full
detail and the cross-condition estimators: `components/laris.md`, `ligand_receptor.md`.

## Step 8 — Downstream: regulons on RNA-only spatial data (cytorete) *(from tutorial)*

```python
import cytorete
piaso.data.fetch_2bit("hg38"); piaso.data.fetch_jaspar(); piaso.data.fetch_genome("hg38")   # cached under ~/.piaso
cytorete.inferRegulon(ds, "hg38", "leiden", jaspar_path=piaso.data.fetch_jaspar(), twobit_path=piaso.data.fetch_2bit("hg38"))
cytorete.pl.regulonEmbedding(ds, regulons=["SOX2", "SPI1"], basis="spatial")                 # activity drawn on the tissue
```
`inferRegulon` is the RNA-only path (promoter motifs × COSG co-specificity) and needs no ATAC;
`inferGRN` (multiome) does not apply to Xenium. See `gene_regulatory_networks.md`.

## Step 9 — Pathways per cell on a section *(from tutorial)*

`matrix, names, _ = piaso.tl.score(adata, gene_list=kegg_sets, layer="infog")` → an AnnData of
cells × pathways with `obsm["spatial"]` copied over → `piaso.pl.embedding(scored, color=[...],
basis="spatial", vmin_pct=2, vmax_pct=98)`. Case-fold gene symbols for mouse libraries (`CUX2` vs
`Cux2`) or you score 1 pathway instead of 277.
