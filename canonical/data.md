# PIASO-data — datasets and reference files for every code block

Every runnable code block in this hub loads from **PIASO-data**, the ecosystem's data repository
(`github.com/genecell/PIASO-data`). It has two halves: **tutorial datasets** on Zenodo and
**genome reference files** committed in the repo. `piaso.data` is the client: it reads the
registry (`datasets.json`, fetched from the repo and cached 24 h), downloads on demand, verifies
md5, and caches under `~/.piaso/data/` (override with `data_dir=`, `piaso.settings.data_dir`, or
`PIASO_DATA_DIR`).

## Zenodo record

- Concept DOI (always the latest version): **10.5281/zenodo.19699638**
- Current record: <https://zenodo.org/records/22012620> (v0.9.0, 2026-08-22; the older record
  19699639 still resolves)
- Direct file URL pattern: `https://zenodo.org/api/records/22012620/files/<filename>/content`

## The registry, from Python

```python
import piaso
piaso.data.list_datasets()                       # table of every dataset with size
info = piaso.data.dataset_info("e18_v3_nuclei")  # dict: title, url, md5, size_bytes, cells, format, counts_layer, tutorials
path = piaso.data.fetch_dataset("e18_v3_nuclei") # download + md5-verify, return the local Path
adata = piaso.data.load_dataset("e18_v3_nuclei") # download + open (AnnData for h5ad / 10x h5, DataFrame for csv)
ds = piaso.data.load_dataset("sea_ad_mtg_20k_cytome", return_type="cytome")  # open cytome Dataset
```

`load_dataset(..., return_type="cytome")` on an h5ad/h5 entry converts once and reopens the file
on later calls. **Check `counts_layer`** before normalizing: `sea_ad_mtg_20k` keeps UMIs in
`layers["UMIs"]`, `adult_cortex_multiome_rna` in `layers["raw"]` (its `.X` is scaled — 93 % of the
values are negative); the 10x `.h5` entries have raw counts in `.X`.

## Datasets (registry v2, 2026-09-04)

| id | format | size | species | cells | counts | reference |
|---|---|---|---|---|---|---|
| `e18_v3_nuclei` | 10x h5 | **19.3 MB** | mouse | ~5,000 | `.X` | 10x Genomics (E18 brain nuclei 5K v3.1) |
| `e18_v3_cell` | 10x h5 | 45.4 MB | mouse | ~10,000 | `.X` | 10x Genomics |
| `e18_v4_cell` | 10x h5 | 64.5 MB | mouse | ~10,000 | `.X` | 10x Genomics |
| `mouse_brain_10k_gemx` | 10x h5 | 65.5 MB | mouse | 11,357 | `.X` | 10x Genomics |
| `pbmc_multiome_san1` | 10x h5 (RNA+ATAC) | 73.2 MB | human | 3,545 | `.X` | De Rop et al. Nat Biotechnol 42, 916–926 (2024) |
| `pbmc_multiome_san2` | 10x h5 (RNA+ATAC) | 83.7 MB | human | 4,360 | `.X` | De Rop et al. (2024) |
| `piaso_markerdb_allen_immune` | csv | **0.1 MB** | human | — | — | Gong et al. Nature 648, 696–706 (2025) |
| `sea_ad_mtg_20k` | h5ad | 1.79 GB | human | 20,000 | `layers["UMIs"]` | Gabitto et al. Nat Neurosci 27, 2366–2383 (2024) |
| `adult_cortex_multiome_rna` | h5ad | 2.48 GB | mouse | 17,412 | `layers["raw"]` | Bravo González-Blas et al. Nat Methods 20, 1355–1367 (2023) |
| `sea_ad_mtg_20k_cytome` | **cytome** | 269 MB | human | 20,000 | `RNA_counts` | Gabitto et al. (2024) |
| `adult_cortex_multiome_rna_cytome` | cytome | 192 MB | mouse | 17,412 | `RNA_counts` | Bravo González-Blas et al. (2023) |
| `allen_devvis_rna` | cytome | 1.41 GB | mouse | 200,061 | `RNA_counts` | Gao et al. Nature 647, 127–142 (2025) |
| `humandevcx_38_rna` | cytome | 1.14 GB | human | 213,090 | `RNA_counts` | Wang et al. Nature 647, 169–178 (2025) |
| `humanlifespan_pfc_rna` | cytome | 25.7 GB | human | 1,501,089 | `RNA_counts` | Catching et al. Cell Reports 45, 117110 (2026) |

Every `.cytome` entry stores raw UMI counts (verified integer at conversion) plus the source
atlas's cell annotations, ready to stream: `piaso.data.load_dataset(name, return_type="cytome")`.

## Fixtures used by the hub's tests

| Purpose | id / file | size | md5 |
|---|---|---|---|
| Loadable AnnData / cytome (real scRNA) | `e18_v3_nuclei` → `SC3_v3_NextGem_DI_Nuclei_5K_..._count_sample_feature_bc_matrix.h5` | 19.3 MB | `81a6ceb41e2def93ac0d0f824a610849` |
| PIASOmarkerDB static slice (offline marker work) | `piaso_markerdb_allen_immune` → `PIASOmarkerDB_AllenHumanImmuneHealthAtlas_L2_251219.csv` | 115 KB | `d4177960c47f995562ad572bb8a5f9f7` |
| **Spatial** (LARIS, image-less) | `adata_tonsil.h5ad` — Slide-tags human tonsil, 5,695 cells × 25,583 genes, `obsm['X_spatial']`, `obs['cell_type']`, from LARIS's Zenodo record **10.5281/zenodo.19981287** (not part of PIASO-data) | 241 MB | — |
| cytorete real run (opt-in, nightly) | `sea_ad_mtg_20k_cytome` + `piaso.data.fetch_jaspar()` + `fetch_2bit("hg38")` (~800 MB) + `fetch_genome("hg38")` | ~1.1 GB | — |

Load the small fixture without `piaso.data` (e.g. in a cosg-only or cytome-only environment):

```bash
curl -L -o e18_v3_nuclei.h5 \
  "https://zenodo.org/api/records/22012620/files/SC3_v3_NextGem_DI_Nuclei_5K_SC3_v3_NextGem_DI_Nuclei_5K_count_sample_feature_bc_matrix.h5/content"
```
```python
import piaso
adata = piaso.pp.read_10x_h5("e18_v3_nuclei.h5")   # raw UMI counts in .X, ~5k nuclei
```

The multi-GB atlases exist for the tutorials that need them (GDR at scale, 1.5M cells) — do
not use them in CI. **PIASO-data still ships no spatial dataset**; the hub's spatial blocks run on
the LARIS tonsil object above, and the Xenium / Stereo-seq / MERFISH blocks are lifted from the
executed piaso.org tutorials (which download from 10x / CNGB / the original providers).

## Genome and motif references (via `piaso.data`)

```python
piaso.data.fetch_genome("hg38")      # gene bodies, promoters, TSS BED, cCREs, chrom sizes (also "mm10"); GTF optional
piaso.data.fetch_2bit("hg38")        # UCSC .2bit genome sequence, ~800 MB, opt-in
piaso.data.fetch_jaspar()            # JASPAR CORE vertebrates MEME
piaso.data.load_lr_database("mouse") # CellChatDB (human 2951 / mouse 3105 pairs) — SCALAR and LARIS input
piaso.data.load_chembl_targets()     # ChEMBL drug-target gene sets
piaso.data.fetch_screen("hg38")      # SCREEN cCRE registry
```

`hg38/` and `mm10/` in the repo hold the ENCODE/UCSC-derived BED files (~17 MB / ~11 MB) that
`fetch_genome` downloads; cytorete needs the TSS BED plus the `.2bit` sequence.

## License

PIASO-data has **no LICENSE file**; its README states genome files derive from public UCSC/ENCODE
annotations and tutorial datasets are **redistributed under CC BY 4.0 with attribution to original
sources**. Cite the original dataset paper (the `reference` column above) when using a fixture.
