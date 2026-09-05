# workflows/

One file per cross-component analysis task. Each opens with its exact multi-package install
line and tracks the AnnData / cytome state passed between steps — this tier is where the
ecosystem is documented as a whole (no single component repo covers composition). Blocks are
executed on the fixtures in `data.md` unless marked *(from tutorial)*.

| File | Task |
|---|---|
| `end_to_end_scrnaseq.md` | load → QC → doublets → INFOG → SVD → Leiden → COSG → artefact check → cell types → GDR (PIASO-native, no scanpy) |
| `streaming_large_data.md` | the same pipeline on a `.cytome`, results written back to the file, R handoff |
| `marker_based_annotation.md` | marker sets (curated / reference) → `predictCellTypeByMarker`; reference projection with `projectGDR`; `predictCellTypeByGDR` |
| `markerdb_annotation.md` | PIASOmarkerDB: `getMarkers` (tuple!) and `analyzeMarkers` — from a gene list to cell types |
| `ligand_receptor.md` | SCALAR (dissociated) vs LARIS (spatial); same CellChatDB |
| `spatial_transcriptomics.md` | Xenium / Visium HD / MERFISH / Stereo-seq → cytome → clusters on the tissue image → ROI → LARIS / cytorete |
| `gene_regulatory_networks.md` | cytorete: regulons from RNA, activity + specificity, plots |
