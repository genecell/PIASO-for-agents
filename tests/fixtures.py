"""Download the smallest PIASO-data fixtures into a local cache for the test suite.

Uses the Zenodo record referenced by genecell/PIASO-data (DOI 10.5281/zenodo.19699639).
Kept tiny on purpose (20 MB h5 + 117 KB csv); the multi-GB atlases are never used in CI.
"""
from __future__ import annotations
import os, urllib.request
from pathlib import Path

CACHE = Path(os.environ.get("PIASO_DATA_CACHE", "/tmp/piaso-data-cache"))

FIXTURES = {
    # id -> (filename on disk, zenodo file name)
    "e18_v3_nuclei": (
        "e18_v3_nuclei.h5",
        "SC3_v3_NextGem_DI_Nuclei_5K_SC3_v3_NextGem_DI_Nuclei_5K_count_sample_feature_bc_matrix.h5",
    ),
    "markerdb_allen_immune": (
        "markerdb_allen_immune.csv",
        "PIASOmarkerDB_AllenHumanImmuneHealthAtlas_L2_251219.csv",
    ),
}
ZENODO = "https://zenodo.org/api/records/19699639/files/{name}/content"


def get(fixture_id: str) -> Path:
    fname, zname = FIXTURES[fixture_id]
    CACHE.mkdir(parents=True, exist_ok=True)
    dest = CACHE / fname
    if not dest.exists() or dest.stat().st_size == 0:
        url = ZENODO.format(name=zname)
        urllib.request.urlretrieve(url, dest)  # noqa: S310 (fixed host)
    return dest


if __name__ == "__main__":
    for k in FIXTURES:
        print(k, "->", get(k))
