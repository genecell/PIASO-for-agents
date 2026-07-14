# AGENTS.md — genecell/PIASO

This repository is part of the **PIASO single-cell omics ecosystem**. Full, cross-component, agent-neutral documentation (with runnable, tested code blocks for every component in Python and R) lives in the hub:
**https://github.com/genecell/PIASO-for-agents**

## Ecosystem at a glance
- **piaso** (`piaso-tools`, python): install `pip install piaso-tools "matplotlib<3.9"`
- **cosg** (`cosg`, python): install `pip install cosg`
- **cosgr** (`library(COSG)`, r): install `remotes::install_github("genecell/COSGR")`
- **laris** (`laris`, python): install `pip install laris`
- **emergene** (`emergene`, python): install `pip install emergene`

## Cross-component decision rules
- **lr_scalar_vs_laris**: Spatial coordinates present (.obsm['spatial']/['X_spatial']; Visium/MERFISH/Xenium) -> LARIS (runLARIS). Dissociated single-cell, no coordinates -> PIASO piaso.tl.runSCALAR. SCALAR needs a user-supplied LR-pair list + specificity matrix (no bundled DB); LARIS bundles CellChatDB.
- **cosg_python_vs_r**: Infer language from objects in session: AnnData / .h5ad / scanpy -> Python `cosg`; Seurat object / .rds / library(Seurat) -> R COSGR. Ask when ambiguous. NOTE default flips: remove_lowly_expressed False(Py)/TRUE(R); n_genes 50(Py)/100(R).
- **emergene_conditions**: >=2 conditions -> runEMERGENE (needs condition_key); single condition -> runMarkG.

For full API, workflows, and citations, read the hub.

---
Maintained by **[The Fishell Laboratory](https://fishelllab.hms.harvard.edu)** (Harvard Medical School / Broad Institute).
