# Changelog

All notable changes to the PIASO-for-agents hub and the `piaso-mcp` server. The hub and the
server are versioned separately: the hub tag is the content release, `piaso-mcp` is the PyPI
package that bundles a snapshot of it.

## Hub v0.2.0 · piaso-mcp 0.1.0 — 2026-09-05

Content refreshed for the August–September 2026 ecosystem releases and **tested against**
piaso-tools 1.2.3 · cosg 1.2.0 · cytome 0.3.1 · laris 0.13.0 · emergene 1.0.2 · cytorete 0.1.1 ·
COSGR 1.0.0 · cytome (R) 0.1.0. The previous release (July 2026) described piaso-tools 1.1.0,
cosg 1.0.4 and laris 0.9.3.

### Highlights

- **Three new components.** [cytome](https://github.com/genecell/cytome) (the single-file
  `.cytome` store every component streams from), **cytome (R)** (read / write / stream the same
  file into Seurat or SingleCellExperiment, no Python), and
  [cytorete](https://github.com/genecell/cytorete) (cell-type-resolved gene regulatory networks /
  regulons on RNA). The ecosystem is documented as four layers: analysis (PIASO), storage
  (cytome), methods (COSG/COSGR, LARIS, Emergene, cytorete), data (PIASO-data).
- **PIASO is documented as a self-contained, scanpy-free pipeline**: `piaso.pp.read_10x_h5` →
  `calculateCellMetrics` → `scrublet` → `filter_cells` → `infog` → `runSVD` → `neighbors` →
  `leiden` → `umap` → `cosg.cosg` → `runGDR`, plus `piaso.data` (datasets, genomes, motif DBs,
  CellChatDB) and the `piaso.pl` suite. scanpy is an optional interop extra.
- **Every function takes an AnnData or a `.cytome`** (`data=`), streams on the latter and writes
  results back onto the file — one contract, documented once, with a naming-contract table.
- **Seven cross-component decision rules**: SCALAR vs LARIS (rewritten — both now take the same
  CellChatDB), COSG Python vs R (defaults updated), Emergene conditions, **AnnData vs cytome**,
  **annotation route** (marker sets / `projectGDR` / `predictCellTypeByGDR` / PIASOmarkerDB),
  **markers vs regulons** (COSG vs cytorete), and **R users** (COSGR + cytome (R) + a `.cytome`
  handoff for everything else).
- **Three new workflows**: streaming on a `.cytome`, spatial transcriptomics (Xenium / Visium HD /
  MERFISH / Stereo-seq → tissue-image overlay → ROI → LARIS → cytorete), gene regulatory networks.
  End-to-end scRNA-seq rewritten PIASO-native with an artefact-cluster check; ligand–receptor and
  both annotation workflows rewritten.
- **The piaso.org tutorial index (46 executed tutorials)** is generated into every target and
  exposed by the MCP server, so an agent can route the user to the executed tutorial for their
  platform before writing code.
- **The `matplotlib<3.9` pin is gone** from every install line; it applied to piaso-tools ≤ 1.1.0
  only. A test now fails if it reappears.

### Behaviour changes an agent will notice

- COSG: `remove_lowly_expressed` now defaults to `True` in Python as in R; `n_genes_user` still
  differs (50 / 100). `cosg.cosg` is polymorphic — on a `.cytome` it **returns a dict** (no
  `key_added`). Analytic p-values (`calculate_pvalues=True`), `batch_key`, GPU documented.
- LARIS 0.13: `runLARIS(lr_data, data, ...)` (`adata=` deprecated), exact p-values via
  `prepareLRBackground` + `background=`, tuple return, new diagnostics; p-values changed in
  0.12 and 0.13.
- SCALAR: the specificity matrix comes from `piaso.tl.specificity_matrix` and the pair list from
  `piaso.data.load_lr_database` (CellChatDB) — the old "user-supplied, no bundled DB" statement is
  gone. The released `runSCALAR` (1.2.3) is the permutation version.
- GDR defaults `n_gene=20`, `mu=10`; `runGDR(save_reference=True)` + `projectGDR` for frozen
  reference spaces; `predictCellTypeByGDR` requires piaso-tools ≥ 1.2.3.
- `runSVDLazy`, `runLARIS(adata=)` and `n_permutations`-only LARIS are documented as deprecated in
  a translation table (`gotchas.md`); no canonical block uses them.
- New gotchas found by execution: `getMarkers(as_dict=True)` returns a tuple; INFOG and COSG
  p-values want integer counts; cytome's `counts` invariant and embedding renames; one writer per
  `.cytome`; `emergene.tl.runEMERGENE` needs `use_rep_acrossDataset`; `leiden_local` on AnnData
  takes `dr_method="X_svd_full"`; `predictCellTypeByGDR` needs disjoint cells and the reference
  label column on the query; `cytorete.inferRegulon` needs `py2bit` (`cytorete[motif]`);
  `cytome.from_10x_h5` needs `h5py`; `regulonSpecificityScatter(key="X_regulon")`.

### piaso-mcp 0.1.0

- New tools: `list_tutorials(topic, component)`, `version_matrix()`, `check_versions()` (current
  PyPI release vs the tested version), `list_datasets()` / `get_dataset(name)` proxying the live
  PIASO-data registry (24 h cache, bundled fallback).
- `resolve_install` knows cytome / cytorete, the R routes (r-universe, conda-forge) and extras;
  `compare_implementations("cytome")` covers the Python vs R cytome packages; `get_api` links the
  generated API reference at https://piaso.org/api/.
- **Compatible with the MCP Python SDK 1.x and 2.x.** SDK 2.0 (2026-07-28) renamed `FastMCP` →
  `MCPServer`; because `mcp>=1.2.0` now resolves to 2.x, **`piaso-mcp 0.0.2` fails on a fresh
  `uvx piaso-mcp`**. 0.1.0 imports either and pins `mcp<3`. Upgrade: `uvx piaso-mcp` picks it up
  automatically; `pip install -U piaso-mcp` otherwise.

### Build, tests, CI

- `build.py` also emits `dist/llms/piaso.org/llms*.txt` (absolute links, drop-in for piaso.org),
  generates `tutorials.md`, stamps every target with "Tested against …", hard-fails if the Claude
  skill description exceeds 1024 characters (1018 now), and no longer truncates the decision rules
  in the Cursor / Copilot files.
- Tests execute every canonical block on the PIASO-native pipeline (37 tests), the cytome
  streaming path, the genome-free cytorete path, and plotting; opt-in heavy tests run LARIS on the
  Slide-tags tonsil (Zenodo 10.5281/zenodo.19981287) and cytorete `inferRegulon` on the SEA-AD
  cortex cytome. Clean-environment checks for cosg (no scanpy), cytome, laris, emergene. R test
  gains an optional cytome (R) round-trip. CI pins the tested versions, runs a 4-way clean-env
  matrix (fixing the nightly failure since 2026-08-28), installs R packages from r-universe, and
  runs the heavy job nightly.

### Known limitations

- Xenium / Stereo-seq / MERFISH blocks are lifted from the executed piaso.org tutorials and marked
  *(from tutorial)*; only the tonsil-based LARIS steps were re-executed here.
- The blind-router activation set was extended (24 prompts) but not re-judged.
- PIASO-data still ships no spatial dataset; the hub uses LARIS's tonsil object for spatial tests.

## Hub v0.1.0 · piaso-mcp 0.0.2 — 2026-07-15

Initial public release: canonical knowledge pack for PIASO 1.1.0, COSG 1.0.4 / COSGR 1.0.0,
LARIS 0.9.3, Emergene 1.0.2 and PIASO-data; generated Claude skill, Cursor rule, Copilot
instructions, AGENTS.md fan-out, llms.txt; `piaso-mcp` local stdio server proxying the live
PIASOmarkerDB API; listed in the official MCP Registry.
