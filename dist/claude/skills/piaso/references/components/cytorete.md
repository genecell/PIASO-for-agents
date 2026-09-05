# cytorete — component reference (self-sufficient)

**cytorete** (*cyto-* + *rete*, "the cell's network"; sy-toh-REE-tee) infers **cell-type-resolved
gene regulatory networks**: transcription factor → target **regulons** from RNA alone, by
combining a **promoter motif scan** (cistrome) with **COSG co-specificity** across cell types,
then scoring each regulon's per-cell **activity** (PIASOscore) and per-type **specificity**
(COSG). It is a SCENIC analogue built on the PIASO stack — COSG replaces GENIE3 co-expression,
the promoter scan replaces RcisTarget pruning, PIASOscore replaces AUCell — and it works on an
AnnData or streams from a `.cytome`. Dependency runs one way: `cytorete → piaso-tools → cosg +
cytome`; cytorete is never a dependency of PIASO. **This release ships the RNA regulon chain**;
the multiome (RNA + ATAC) GRN chain and the ATAC TF-activity chain exist as names that raise
`ImportError` at call time. Tutorials: https://piaso.org/tutorials/cytorete-regulons/,
`cytorete-spatial-stereoseq`, `cytorete-regulon-dynamics`, `motif-analysis`.

## Install

```bash
pip install cytorete            # pulls piaso-tools>=1.2.2, cosg, cytome
pip install 'cytorete[motif]'   # + py2bit — REQUIRED to read genome sequence from a .2bit, i.e. for inferRegulon
```

Executed against **cytorete 0.1.1** with piaso-tools 1.2.3. Without `py2bit`, `inferRegulon` stops
with `ImportError: piaso.data needs the optional 'py2bit' package ...` after the COSG step —
install the extra up front. BSD-3-Clause; no paper yet — cite the repository.

## Import

```python
import cytorete                 # cytorete.tl (== .tools), .pp, .pl, .data; camelCase and snake_case both exported
cytorete.inferRegulon           # == cytorete.tl.inferRegulon == cytorete.infer_regulon
```

The `piaso.tl.inferRegulon` / `regulonActivity` / `regulonSpecificity` / `inferGRN` /
`inferTFActivity` names are **thin forwarders** to cytorete (they raise `ImportError` pointing at
`pip install cytorete` if it is missing) — new code should import cytorete directly.

## Reference data (all fetched and cached by `piaso.data`)

```python
import piaso
jaspar = piaso.data.fetch_jaspar()          # JASPAR 2024 CORE vertebrates, MEME format
twobit = piaso.data.fetch_2bit("hg38")      # UCSC genome sequence, ~800 MB, opt-in  (also "mm10")
piaso.data.fetch_genome("hg38")             # TSS BED (+ promoters, cCREs) among the PIASO-data references
```

`motif_db="cisbp"` / `"both"` with `piaso.data.fetch_cisbp()`; `tf_list_path=` or
`piaso.data.fetch_animaltfdb_tf_list("human")` to restrict TFs; `regulatory_regions="promoter+cre"`
adds SCREEN cCREs (`piaso.data.fetch_screen`).

## What `inferRegulon` computes

1. **COSG λ specificity matrix** over all genes for `groupby` (`piaso.tl.specificity_matrix`;
   streamed on a cytome) → candidate target genes (`target_genes='cosg'`).
2. **Promoter cistrome** — strand-aware windows (`upstream=1000`, `downstream=500` around each
   TSS; alternative promoters kept separate; `biotypes=('protein_coding',)`) extracted from the
   `.2bit`, scanned against the motif DB with `piaso.pp.scan_motifs` (Rust) and a per-motif
   background (`cistrome_method='motif_bg'`; `'nes'` gives RcisTarget-style NES pruning) → TF × gene
   motif support.
3. **Trans co-specificity** (`piaso.tl.cospecificity_trans`) — each motif-supported TF→gene pair is
   kept when the TF's and target's COSG profiles agree (cosine, positive sign) across cell types.
4. **Regulon assembly** — TFs with ≥ `min_targets` surviving targets; global and per-cell-type regulons.
5. **Activity** — `piaso.tl.score` per regulon per cell against expression-matched control sets,
   with per-cell p-values (`compute_pvalues=True` default, `n_ctrl_set=1000`).
6. **Specificity** — COSG on the activity matrix → regulon × cell type (RSS analogue).

**Reads** raw counts (normalizes with INFOG internally; `cosg_layer` / `score_layer` if you have
layers), `obs[groupby]` / `cells[groupby]`. **Writes** onto the object: embeddings **`X_regulon`**
(cells × regulons, activity) and **`X_regulon_pval`**, and a `regulon` entry in `uns` /
`ds.metadata` with keys `regulons`, `weights`, `per_celltype`, `edges`, `cistrome_density`,
`celltypes`, `params`, `names`, `activity_key`, `tf_pct`, `specificity`, `specificity_matrix`.
Returns `None` unless `copy=True`.

## Executed block (SEA-AD cortex cytome, 20,000 nuclei, 15 TFs)

```python
import numpy as np, piaso, cytorete
ds = piaso.data.load_dataset("sea_ad_mtg_20k_cytome", return_type="cytome")   # human cortex, 'Subclass' annotations, streamed
TFS = ["SOX9", "OLIG2", "SPI1", "IRF8", "MEF2C", "NEUROD2", "DLX2", "LHX6",
       "TCF7L2", "NFIB", "RORB", "CUX2", "FEZF2", "PAX6", "EGR1"]
cytorete.inferRegulon(ds, "hg38", "Subclass", jaspar_path=jaspar, twobit_path=twobit, tf_list=TFS)   # omit tf_list to run every expressed TF with a motif
# log: 2552 genes with promoters (5393 intervals); cistrome 15 TFs x 2552 genes, 3411 edges (8.9 %);
#      15 global regulons (median 45 targets); per-cell-type regulons for 23 cell types; ~50 s after loading
md = ds.metadata["regulon"]; names = md["names"]          # ['CUX2', 'DLX2', 'EGR1', 'FEZF2', 'IRF8', 'LHX6', ...]
A = np.asarray(ds.embeddings["X_regulon"]); P = np.asarray(ds.embeddings["X_regulon_pval"])   # (20000, 15) each
spec = cytorete.regulonSpecificity(ds, groupby="Subclass", copy=True)     # LONG-form DataFrame: cell_type, regulon, cosg_score, rank — pivot for a heatmap
```

On an AnnData the same call writes `adata.obsm["X_regulon"]`, `adata.obsm["X_regulon_pval"]`,
`adata.uns["regulon"]`. `regulonActivity(data, regulons={TF: [targets]}, weights=None,
score_layer="infog", compute_pvalues=True, key_added="X_regulon")` scores any regulon dict you
supply — including one you did not infer — via `piaso.tl.score`.

## Reading the result — `cytorete.pl`

```python
cytorete.pl.regulonActivity(ds, groupby="Subclass")                       # regulon x cell type heatmap (activity or specificity; style=, values='zscore', significance=)
cytorete.pl.regulonEmbedding(ds, regulons=["CUX2", "FEZF2", "LHX6", "IRF8"], basis="X_umap")   # where is a regulon active (use_pval=True for -log10 p)
cytorete.pl.regulonNetwork(ds, tf="CUX2", max_targets=20)                 # the TF and its targets (label_targets=False for dense layouts)
cytorete.pl.regulonSpecificityScatter(ds, groupby="Subclass", key="X_regulon")   # specificity vs the TF's own expression (its default key is X_grn — pass key=)
```

The controls fall where they should on the cortex: SPI1/IRF8 in microglia, OLIG2 in OPC /
oligodendrocytes, SOX9/PAX6 in astrocytes, LHX6 in MGE-derived interneurons, FEZF2 in deep-layer
projection neurons, CUX2 in L2/3 IT. **Read activity and p-value together**: with 1,000 control
sets −log10 p saturates at 3.0 ("as significant as the test can report"); a regulon confined to a
small population (IRF8, microglia) is highly specific yet significant in few cells — specificity
and per-cell significance are different questions. Sanity-check a regulon by its targets
(CUX2 → GRIA3, KCNQ5, DLGAP1, CDH9) before believing it.

## Spatial and developmental data

Same call on a spatial cytome (`cytorete.inferRegulon("section.cytome", "mm10", "annotation",
jaspar_path=..., twobit_path=piaso.data.fetch_2bit("mm10"), tf_list=..., cosg_layer="infog",
score_layer="infog")`), then `cytorete.pl.regulonEmbedding(ds, regulons=[...], basis="spatial")`
draws regulon activity on the tissue; the Stereo-seq tutorial compares against published regulons
(Spearman per TF). Regulons need the TF **and** its targets measured — check panel coverage first
(`workflows/spatial_transcriptomics.md`). `inferGRN` (multiome) does not apply to RNA-only data.

## Decision rule — markers vs regulons

"Which genes define this type" → **COSG**. "Which TFs drive it" → **cytorete**. `inferRegulon`
needs a genome, a TSS annotation and a motif DB (`hg38` / `mm10`), a cell-type column, and
`py2bit`; it does not need ATAC.

## Withheld chains (this release)

`cytorete.inferGRN`, `inferGRN_consensus`, `inferTFActivity`, `pp.build_peak_cistrome` exist in
the package and raise `ImportError: ... not part of this distribution` at call time, so
`import cytorete` behaves the same either way. `pp.build_cistrome(promoter_seqs, tf_motif_map)`
and `pp.extract_promoter_sequences` are the public pieces of the RNA chain.

## Citation

No paper yet. Cite the repository: https://github.com/genecell/cytorete (Min Dai, The Fishell
Laboratory; BSD-3-Clause). Its dependencies carry their own citations (COSG — Dai et al. 2022;
PIASO — Wu, Dai et al. Nature 2026).
