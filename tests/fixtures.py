"""Download the fixtures the test suite needs into a local cache.

Small fixtures come from the PIASO-data Zenodo record (concept DOI 10.5281/zenodo.19699638,
current record 22012620; the older record 19699639 still resolves). Kept tiny on purpose
(19 MB h5 + 115 KB csv); the multi-GB atlases are never used in CI.

Heavy fixtures (only when PIASO_HEAVY=1): the LARIS Slide-tags tonsil object (241 MB, Zenodo
10.5281/zenodo.19981287) for the spatial tests, and the SEA-AD cortex cytome + JASPAR + hg38 .2bit
(via piaso.data) for the real cytorete run.
"""
from __future__ import annotations
import os, urllib.request
from pathlib import Path

CACHE = Path(os.environ.get("PIASO_DATA_CACHE", "/tmp/piaso-data-cache"))
HEAVY = os.environ.get("PIASO_HEAVY") == "1"

FIXTURES = {
    # id -> (filename on disk, zenodo record, zenodo file name)
    "e18_v3_nuclei": (
        "e18_v3_nuclei.h5", "22012620",
        "SC3_v3_NextGem_DI_Nuclei_5K_SC3_v3_NextGem_DI_Nuclei_5K_count_sample_feature_bc_matrix.h5",
    ),
    "markerdb_allen_immune": (
        "markerdb_allen_immune.csv", "22012620",
        "PIASOmarkerDB_AllenHumanImmuneHealthAtlas_L2_251219.csv",
    ),
    # heavy — LARIS tonsil (spatial), not part of PIASO-data
    "laris_tonsil": ("adata_tonsil.h5ad", "19981287", "adata_tonsil.h5ad"),
}
ZENODO = "https://zenodo.org/api/records/{record}/files/{name}/content"


def get(fixture_id: str) -> Path:
    fname, record, zname = FIXTURES[fixture_id]
    CACHE.mkdir(parents=True, exist_ok=True)
    dest = CACHE / fname
    if not dest.exists() or dest.stat().st_size == 0:
        url = ZENODO.format(record=record, name=zname)
        urllib.request.urlretrieve(url, dest)  # noqa: S310 (fixed host)
    return dest


if __name__ == "__main__":
    for k in FIXTURES:
        if k == "laris_tonsil" and not HEAVY:
            continue
        print(k, "->", get(k))
