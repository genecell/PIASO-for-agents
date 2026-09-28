# PIASO for Claude

A skill that teaches Claude the PIASO single-cell omics ecosystem, so it writes correct analysis code for scRNA-seq and spatial data in Python and R: quality control, INFOG normalization, clustering, marker genes with COSG and COSGR, gene-set scoring, cell-type annotation with PIASOmarkerDB, single-cell (SCALAR) and spatial (LARIS) ligand-receptor analysis, gene regulatory networks with cytorete, differential expression with Emergene, and out-of-core analysis of large datasets stored as `.cytome` files.

## What it contains

Only text: `skills/piaso/SKILL.md`, which routes a request to the right component, and reference pages under `skills/piaso/references/` for each component and cross-component workflow. The plugin runs no code, starts no server, and sends nothing anywhere. The analysis code Claude writes runs in your own environment, with the packages you install (`pip install piaso-tools`).

## Links

- Documentation and executed tutorials: https://piaso.org/tutorials/
- API reference: https://piaso.org/api/
- Source of this plugin: https://github.com/genecell/PIASO-for-agents

Maintained by The Fishell Laboratory (Harvard Medical School / Broad Institute). License: BSD-3-Clause.
