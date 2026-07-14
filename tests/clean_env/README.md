# Clean-env single-component proof

Each script runs one component's canonical workflow in an environment where **only that
component is installed** (no `piaso-tools`), proving `components/<name>.md` is self-sufficient.
Each asserts `import piaso` fails, then runs the documented block on a PIASO-data fixture.

    python -m venv /tmp/cosgonly && /tmp/cosgonly/bin/pip install cosg igraph leidenalg "matplotlib<3.9"
    /tmp/cosgonly/bin/python tests/clean_env/run_cosg.py /path/to/e18_v3_nuclei.h5

(See .github/workflows/test.yml → clean-env-selfsufficiency for the CI version.)
