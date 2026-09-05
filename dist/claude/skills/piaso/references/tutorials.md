# piaso.org tutorials — index

Executed, human-reviewed tutorials at https://piaso.org/tutorials/ (generated API reference: https://piaso.org/api/). Each runs against the real dataset it names; the numbers and figures are what the code produced. Grouped by topic; `components` says which packages it uses.

## scRNA-seq

- [Human PBMC scRNA-seq end to end (AnnData)](https://piaso.org/tutorials/san2-end-to-end-anndata/) — counts -> two-tailed QC -> doublets -> clusters -> artefact clusters -> cell types; the place to start *(piaso, cosg)*
- [Human PBMC scRNA-seq end to end (cytome)](https://piaso.org/tutorials/san2-end-to-end-cytome/) — the same analysis streamed from disk; filter writes a new file; annotation written back to the file *(piaso, cosg, cytome)*
- [Mouse brain scRNA-seq end to end (AnnData)](https://piaso.org/tutorials/rna-end-to-end-anndata/) — Cell Ranger counts -> QC -> doublets -> clusters -> two annotation routes (PIASOmarkerDB vs reference) -> GDR *(piaso, cosg)*
- [Mouse brain scRNA-seq end to end (cytome)](https://piaso.org/tutorials/rna-end-to-end-cytome/) — the same mouse analysis, streamed *(piaso, cosg, cytome)*
- [Multiple samples: QC, embedding and annotation](https://piaso.org/tutorials/rna-multi-sample/) — per-library doublets (scrublet library_key), per-sample QC, GDR with batch_key vs SVD, silhouette check, runHarmony when a batch effect is real *(piaso, cosg)*
- [Two preparations of one tissue: cells vs nuclei](https://piaso.org/tutorials/rna-multi-sample-mouse/) — a real batch effect, measured, and where it comes from *(piaso, cosg)*
- [Multiple samples in one cytome](https://piaso.org/tutorials/rna-multi-sample-cytome/) — cytome.merge, then the multi-sample analysis streamed from one file *(piaso, cytome)*

## methods

- [Marker-gene-guided dimensionality reduction (GDR)](https://piaso.org/tutorials/gdr/) — runGDR / runGDRParallel with batch_key and groupby=None; batch mixing vs cell-type purity measured *(piaso)*
- [Cell type prediction with GDR](https://piaso.org/tutorials/gdr-predict-cell-type/) — predictCellTypeByGDR: label transfer from an annotated reference *(piaso)*
- [projectGDR: putting new data into a reference's space](https://piaso.org/tutorials/project-gdr/) — runGDR(save_reference=True) then projectGDR(query, reference); mode='reference' vs 'self'; novelty flags *(piaso)*
- [Local Leiden clustering](https://piaso.org/tutorials/leiden-local/) — piaso.tl.leiden_local: re-cluster inside selected groups *(piaso)*
- [GDR at scale: 200,000 cells in 17 minutes](https://piaso.org/tutorials/gdr-at-scale/) — the whole pipeline on a .cytome path; memory/time knobs *(piaso, cytome)*
- [GDR and SVD on 1.5 million cells](https://piaso.org/tutorials/gdr-vs-svd-1p5m/) — runtime, peak memory and separation on a 1.5M-nucleus atlas *(piaso, cytome)*
- [GDR on developmental data](https://piaso.org/tutorials/gdr-developmental/) — stages rather than terminal types *(piaso)*
- [GDR beyond transcriptomics](https://piaso.org/tutorials/gdr-applications/) — images as expression matrices *(piaso)*
- [Emergene: differential analysis at the level of the cell](https://piaso.org/tutorials/emergene/) — runEMERGENE downstream of PIASO (use_rep='X_gdr', layer='infog'); per-condition vs all-in-one *(emergene, piaso)*

## marker-genes

- [COSG: marker genes and their significance](https://piaso.org/tutorials/cosg-markers/) — mu / expressed_pct / n_genes_user; calculate_pvalues columns; plotMarkerDotplot / plotMarkerDendrogram; IQR normalisation across cell types; the double-dipping caveat *(cosg, piaso)*
- [COSG on a cytome](https://piaso.org/tutorials/cosg-cytome/) — run_cosg_cytome / cosg.cosg(ds); output_format shapes; what layer= means on a file *(cosg, cytome)*
- [COSG across batches](https://piaso.org/tutorials/cosg-batch/) — batch_key scores per batch and averages; batch_cell_number_threshold *(cosg)*
- [COSG on the GPU](https://piaso.org/tutorials/cosg-gpu/) — device='gpu'; measured speed-ups by matrix size *(cosg)*
- [COSG on spatial transcriptomics](https://piaso.org/tutorials/cosg-spatial/) — organ markers on a whole embryo section, plotted in tissue space *(cosg, piaso)*

## annotation

- [Marker-based cell type prediction](https://piaso.org/tutorials/marker-cell-type-prediction/) — predictCellTypeByMarker from COSG markers (mu=10) on a GDR embedding; 95.4% held-out accuracy; the same call against PIASOmarkerDB *(piaso, cosg)*
- [PIASOmarkerDB API client](https://piaso.org/tutorials/markerdb-api/) — getMarkers by study/gene/cell_type/species; as_dict=True returns a tuple; analyzeMarkers from a gene list to cell types *(piaso)*

## gene-sets

- [Gene set scoring (PIASOscore)](https://piaso.org/tutorials/gene-set-scoring/) — piaso.tl.score with compute_pvalues; why matched controls; a whole pathway database at once; COSG on the score matrix; on a cytome *(piaso, cosg)*
- [KEGG and drug-target gene sets](https://piaso.org/tutorials/kegg-drug-targets/) — 320 pathways + 659 drug target sets scored per cell, then which cell type each belongs to *(piaso, cosg)*
- [Motif analysis](https://piaso.org/tutorials/motif-analysis/) — genome -> promoters -> PWM scan (piaso.pp.scan_motifs) -> enrichment; the background choice *(piaso)*

## spatial

- [Xenium with tissue-image overlay](https://piaso.org/tutorials/spatial-xenium/) — read_10x_h5 + cells.csv -> obsm['spatial'] -> infog/SVD/leiden -> cytome.from_anndata + add_spatial_image -> plotEmbedding(image=True) -> cells_in_region ROI *(piaso, cytome, cosg)*
- [Xenium into a cytome directly](https://piaso.org/tutorials/spatial-xenium-cytome/) — cytome.from_10x_h5 -> add_embedding('spatial') + set_spatial_coords -> streamed clustering -> image overlay *(cytome, piaso)*
- [Xenium Prime 5K: mouse brain end to end](https://piaso.org/tutorials/xenium-prime-brain/) — 63,173 cells x 5,006 genes from the raw bundle to annotated clusters *(piaso, cosg, cytome)*
- [Atera WTA: near-whole-transcriptome in situ](https://piaso.org/tutorials/atera-wta/) — 18,028 targets in situ *(piaso, cosg, cytome)*
- [Downstream in situ: pathways, ligand-receptor and regulons](https://piaso.org/tutorials/xenium-downstream/) — LR-pair coverage check first; KEGG per cell; LARIS on a section (plotCCCSpatial); cytorete.inferRegulon on RNA-only spatial data *(piaso, laris, cytorete)*
- [MERFISH sections in one cytome](https://piaso.org/tutorials/spatial-merfish/) — several sections in one file, plot_embeddings_split by section *(cytome, piaso)*
- [Stereo-seq whole embryo](https://piaso.org/tutorials/spatial-stereoseq/) — bin-resolution section; rotateSpatialCoordinates on a cytome; zoom by coordinates (cell_mask) *(cytome, piaso)*
- [GDR on spatial transcriptomics: eight embryonic stages](https://piaso.org/tutorials/gdr-mosta-allstage/) — 520,815 bins in one embedding *(piaso, cytome)*

## grn

- [cytorete: RNA regulon inference](https://piaso.org/tutorials/cytorete-regulons/) — fetch_jaspar/fetch_2bit/fetch_genome -> inferRegulon(ds, 'hg38', 'Subclass', tf_list=) -> X_regulon(+_pval) -> regulonActivity / regulonSpecificity / four plots *(cytorete, piaso, cytome)*
- [cytorete on spatial data](https://piaso.org/tutorials/cytorete-spatial-stereoseq/) — regulons across a whole embryo section, compared against published regulons *(cytorete, cytome, piaso)*
- [Regulon dynamics across development](https://piaso.org/tutorials/cytorete-regulon-dynamics/) — half a million bins, eight stages *(cytorete, cytome)*

## cell-cell-interaction

- [SCALAR: ligand-receptor interaction analysis](https://piaso.org/tutorials/scalar/) — load_lr_database -> COSG all-gene specificity matrix -> runSCALAR (exact matched-gene null, 1.67M interactions) -> plotLigandReceptorInteraction / Lollipop; split by CellChatDB annotation *(piaso, cosg)*
- [LARIS: ligand-receptor interactions in spatial data](https://piaso.org/tutorials/laris/) — SCALAR vs LARIS decision; links the LARIS repo tutorials 01-08 *(laris)*

## plotting-data

- [Plotting](https://piaso.org/tutorials/plotting/) — piaso.pl end to end: embedding, dotplot, violin, scatter, split embeddings *(piaso)*
- [Colour palettes](https://piaso.org/tutorials/color-palettes/) — built-in palettes (piaso.pl.color.*) and custom ones *(piaso)*
- [Datasets and genome references](https://piaso.org/tutorials/datasets-and-references/) — piaso.data.list_datasets / load_dataset / fetch_dataset / dataset_info; genomes; caching and data_dir *(piaso)*
- [cytome basics](https://piaso.org/tutorials/cytome-basics/) — create / from_anndata / from_h5ad(backed=True) / from_cellranger; cells.query; RNA.counts slicing and iter_rows; embeddings; provenance; set_categories *(cytome, piaso)*
- [Converting: AnnData, Seurat, SingleCellExperiment](https://piaso.org/tutorials/cytome-conversion/) — the six conversions; what travels (graphs by default, layers opt-in, feature annotations beyond id/symbol do not) *(cytome, cytome-r)*
- [cytome in R](https://piaso.org/tutorials/cytome-in-r/) — read_cytome(as=), write_cytome, delayed=TRUE out-of-core, cytome_stream *(cytome-r)*
- [Agents and project tooling](https://piaso.org/tutorials/agents-and-tooling/) — PIASO-for-agents, stato, PlanDrop — why an agent should read the current API *(—)*
