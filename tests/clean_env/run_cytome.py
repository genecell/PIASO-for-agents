"""cytome-only environment (pip install cytome h5py): no piaso, no cosg, no scanpy — build a .cytome
from the 10x h5, query it, stream it, write a column. Proves components/cytome.md is self-sufficient.
(h5py is needed by cytome.from_10x_h5 but not declared by cytome 0.3.1 — see gotchas.md.)"""
import warnings; warnings.filterwarnings("ignore")
import sys, os, tempfile
for mod in ("piaso", "cosg", "scanpy"):
    try:
        __import__(mod); print(f"FAIL: {mod} present in cytome-only env"); sys.exit(1)
    except ImportError:
        pass
import numpy as np, cytome

path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/piaso-data-cache/e18_v3_nuclei.h5"
out = os.path.join(tempfile.mkdtemp(), "e18.cytome")
ds = cytome.from_10x_h5(path, out, sample_name="E18")           # use the returned (open) dataset
assert ds.n_cells > 5000 and "RNA_counts" in ds.list_matrices()
n = 0
for start, end, chunk in ds.RNA.counts.iter_rows():
    n += chunk.shape[0]
assert n == ds.n_cells
depth = np.asarray(ds.RNA.counts[:, :].sum(axis=1)).ravel() if ds.n_cells < 10000 else None
ds.cells["depth"] = depth; ds.flush()
mask = ds.cells.query_mask("depth >= 1000")
assert mask.sum() > 0
ds.close()
ds2 = cytome.open(out); assert "depth" in ds2.cells.columns; ds2.close()
print("cytome-only OK:", n, "cells streamed, column written and read back (no piaso/cosg/scanpy)")
