# LARIS — component reference (self-sufficient)

LARIS (**L**igand **A**nd **R**eceptor **I**nteraction in **S**patial transcriptomics) infers
ligand–receptor interactions for **spatial** data (Slide-tags, Visium / Visium HD, Xenium, MERFISH,
Stereo-seq, …). Interaction strength is computed **per cell against its spatial neighbours**, so a
pair only scores where the partners are adjacent; from that it derives spatially specific LR
pairs, **sender → receiver cell-type scores** with **exact p-values** against an expression-matched
null, spatial-neighbourhood context, and (since 0.10) validated **cross-condition comparison**. It
keys off a spatial kNN graph and is meaningless without coordinates. Standalone package (pulls
`cosg` and scanpy; PIASO is not required, only its figure style is used if present). Executed
against **laris 0.13.0** on the Slide-tags human tonsil (`adata_tonsil.h5ad`, 5,695 cells ×
25,583 genes, Zenodo 10.5281/zenodo.19981287). Eight executed tutorials live in the repo:
https://github.com/genecell/LARIS/blob/master/tutorials/README.md.

## Install

```bash
pip install laris            # compiled Rust kernels in the wheel since 0.12.0; NumPy fallback for source builds (LARIS_NO_RUST=1)
pip install 'laris[cytome]'  # streaming from a .cytome (already satisfied if piaso-tools is installed)
```

## Import / entry points

`import laris as la`: `la.tl` (`prepareLRInteraction`, `prepareLRBackground`, `runLARIS`,
`compareLARIS`, `compareLARISMatched`, `buildJointEmbedding`, `combineComparisons`,
`permuteLRPairs`, `computeDecoyFDR`, `decoyReport`, `readCytome`), `la.pl` (`plotCCCHeatmap`,
`plotCCCDotPlot`, `plotCCCNetwork`, `plotCCCNetworkCumulative`, `plotCCCSpatial`, `plotLRDotPlot`,
`plotCompareLARIS`), `la.datasets.lrDatabase`. Every `data`/`lr_data` argument accepts an AnnData,
an open cytome Dataset or a `.cytome` path (`adata=` / `lr_adata=` are deprecated aliases).

## Bundled database

`la.datasets.lrDatabase(species="human" | "mouse", pathway=None, annotation=None)` → **CellChatDB**
(human 2951 / mouse 3105 pairs; columns `interaction_name`, `ligand`, `receptor`, `pathway_name`,
`annotation` ∈ {Secreted Signaling, ECM-Receptor, Cell-Cell Contact, Non-protein Signaling}, …).
The same tables `piaso.data.load_lr_database` fetches for SCALAR — a pair list curated for one
transfers to the other. Any DataFrame with `ligand` / `receptor` columns works.

## What it computes

1. **`prepareLRInteraction(data, lr_df, use_rep_spatial="X_spatial", number_nearest_neighbors=20,
   sigma="adaptive", unmatched="drop")`** — diffuses ligand and receptor expression over the spatial
   kNN graph (adaptive kernel bandwidth, so coordinate units do not matter) and multiplies the
   diffused pair. **Returns a new AnnData** (cells × LR pairs; `var_names = "ligand::receptor"`;
   `obsm` incl. `X_spatial` carried over). Pairs whose genes are absent are dropped with a warning
   (`unmatched="error"` for a custom database). Requires coordinates in **`obsm["X_spatial"]`**
   (LARIS's key, not scanpy's `spatial`).
2. **`prepareLRBackground(data, lr_df, use_rep_spatial=, n_matched_genes=100, n_pool=4000,
   augment_pool=True)`** — the null: for each LR gene, expression-matched (mean/variance) partner
   genes drawn from a candidate pool that is grown until every LR gene has matched peers above it;
   every pseudo-pair is scored through the full pipeline. Depends only on the cells, graph and
   gene set → **build once, reuse** across every `groupby` and parameter sweep (`pickle` it).
   This is the expensive step (≈ 6 min on the 5,695-cell tonsil with `n_matched_genes=30`; the run below).
3. **`runLARIS(lr_data, data, use_rep="X_spatial", use_rep_spatial="X_spatial", groupby="CellTypes",
   background=bg, n_top_lr=4000, mu=0.25, by_celltype=True, min_null_support=0)`** — step 1: spatial
   specificity of each LR pair (closed-form; deterministic, no seed). Step 2 (`by_celltype=True`):
   COSG specificity of ligand / receptor per cell type × spatial co-localization → per
   sender → receiver interaction scores; with `background=` the p-value is the **exact tail**
   `(#pseudo-pairs ≥ observed + 1) / (n + 1)` (floor 1/(n_matched_genes² + 1); 1e-4 at defaults),
   BH FDR per group. **Returns a tuple** `(laris_lr, res)`.

Output `res` columns: `sender, receiver, interaction_score, ligand, receptor, interaction_name,
p_value, p_value_fdr, nlog10_p_value_fdr` + diagnostics `null_matchability` (≥ 0.99 ⇒ the gene
sits above its whole matched set — p overstated, a warning names them), `null_support`,
`pair_breadth` (fraction of the sender–receiver grid a pair is called in; > ~0.25 ⇒ ubiquity, not
cell-type information).

## Executed block (tonsil)

```python
import anndata as ad, laris as la
adata = ad.read_h5ad("adata_tonsil.h5ad")                       # obsm['X_spatial'], obs['cell_type'], log-normalized .X
lr_df = la.datasets.lrDatabase(species="human")                   # 2951 pairs
lr_data = la.tl.prepareLRInteraction(adata, lr_df, use_rep_spatial="X_spatial")   # (5695, 1985): 1985 pairs have both genes measured
bg = la.tl.prepareLRBackground(adata, lr_df, use_rep_spatial="X_spatial")         # ~5 min; the reusable part
laris_lr, res = la.tl.runLARIS(lr_data, adata, use_rep="X_spatial", use_rep_spatial="X_spatial",
                               groupby="cell_type", background=bg)
res[res.p_value_fdr < 0.05].nsmallest(5, "p_value")[["sender", "receiver", "interaction_name", "interaction_score", "p_value_fdr"]]
```

Executed here: 389,060 sender–receiver–pair combinations, **1,372 significant at FDR < 0.05** with
`n_matched_genes=30` (the tutorial's default `n_matched_genes=100` run reports 1,345); `laris_lr` is a
(1985 × 4) table `ligand, receptor, score, Rank`. With 0.13.0 defaults a plain call reproduces the
published tonsil reference values (`FCER2::CR2` B_naive → B_naive on top by score; `C3::CR2` MRC → FDC
and `SEMA7A::PLXNC1` B_memory → mDC among the most significant). Group by **clusters or coarse types**: LARIS tests every sender–receiver
combination (223 fine types ≈ 50,000 pairs vs 34 clusters ≈ 1,156; hours vs ninety seconds).

### Plots

```python
la.pl.plotCCCHeatmap(res)                                                          # sender x receiver counts of significant pairs (filter_significant=True, threshold=0.05); no show= kwarg on la.pl functions — use plt.close()
la.pl.plotCCCDotPlot(res, interactions_to_plot=["C3::CR2", "FCER2::CR2"], senders=["MRC", "B_naive"], receivers=["FDC_LZDZ", "B_naive"])   # senders/receivers (or sender_receiver_pairs=) are REQUIRED
la.pl.plotCCCNetwork(res, cell_type_of_interest="FDC_LZDZ", interaction_direction="sending", data=adata, groupby="cell_type")   # self-interactions drawn as rings (0.12)
la.pl.plotCCCSpatial(lr_data, basis="X_spatial", interaction="C3::CR2", color_by="score")    # the per-cell score ON the section — what LARIS is for
la.pl.plotLRDotPlot(la.tl.prepareDotPlotAdata(...), interactions_to_plot=[...], groupby="cell_type")
```

### Cross-condition comparison (0.10+, from the LARIS tutorials 02/03/06)

`la.tl.compareLARIS(results, conditionMap=, referenceCondition=, sampleToSubject=, level="both")`
compares per-sample `res` tables with **subject-level** statistics (aggregate estimator; volcano via
`la.pl.plotCompareLARIS`); `compareLARISMatched` + `buildJointEmbedding(adata, batch_key=,
method="harmony" | "gdr" | "pca")` compares **at matched cell states** so a difference is not a
difference in composition; `combineComparisons` merges; `section_key=` handles several sections in
one object; crossed labels (cell type × region) via PIASO's `getCrossCategories`. Decoy control:
`permuteLRPairs` + `decoyReport` measure what the database contributes over chance pairing.

### On a `.cytome`

`prepareLRInteraction("section.cytome", lr_df, output="lr.cytome")`, `runLARIS(lr_ds, "section.cytome",
...)` and `prepareLRBackground` all stream; LARIS tutorial 04 shows the on-disk path is
bit-identical to the in-memory one.

## Version notes an agent will meet

- **0.10**: defaults `mu 1 → 0.25`, `sigma 100 → 'adaptive'`, neighbourhood sizes `10 → 20`,
  `spatial_weight 1 → 3.0` — scores differ from < 0.10 by design.
- **0.12**: spatial specificity became analytic (deterministic; `n_repeats` ignored); exact
  matched-gene null introduced (`prepareLRBackground`); Rust kernels in the wheel.
- **0.13**: the matched pool now covers the whole expression range — **p-values for highly
  expressed genes changed**; some 0.12 calls are no longer significant. `n_permutations` is
  legacy-only (raises `FutureWarning` with `background=`). If you report p-values, read tutorial 07.
- Deprecated aliases still accepted: `adata=` / `lr_adata=`, `run..(n_permutations=)` without a
  background (the old sampled null, kept one cycle).

## Decision rule — LARIS vs SCALAR

Both answer "which ligand–receptor interactions / cell–cell communications are happening?":
**spatial coordinates present** → LARIS; **dissociated single-cell** → SCALAR (`piaso.tl.runSCALAR`,
`components/piaso.md`). Same CellChatDB either way; SCALAR gives one score per (pair, sender,
receiver), LARIS additionally a score per cell that can be drawn on the tissue. Natural pairing:
SCALAR on the dissociated reference, LARIS on the section. On a targeted panel check how many
pairs have both genes measured before either.

## Citation

LARIS is a **preprint** (bioRxiv, 2025) — a separate publication from the *Nature* (2026) paper.

> M. Dai, T. Török, D. Sun, et al. LARIS enables accurate and efficient ligand and receptor
> interaction analysis in spatial transcriptomics. *bioRxiv* (2025). DOI: 10.1101/2025.11.26.690796
