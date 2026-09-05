"""piaso-mcp — a local stdio MCP server for the PIASO ecosystem.

Serves KNOWLEDGE and PUBLIC DATA only — never compute, never user data.
- Doc tools (search_docs, get_api, compare_implementations, resolve_install, list_tutorials,
  version_matrix) read a bundled snapshot of canonical/ (populated by build.py).
- Registry tools (list_datasets, get_dataset) PROXY the live PIASO-data registry
  (datasets.json on GitHub) with the bundled data.md as the offline fallback.
- check_versions asks PyPI for the current release of each package and compares it with the
  versions this snapshot was tested against, so an agent can disclose drift.
- Marker tools (query_marker_db, get_markers, list_studies) PROXY the live PIASOmarkerDB REST API
  at https://piaso.org/piasomarkerdb/api/v1/ (decision: proxy the existing API rather than bundle
  the data). This gives R / non-Python callers their first programmatic path to PIASOmarkerDB.

Needs NO PIASO packages installed. Stdlib + the `mcp` SDK (1.x or 2.x) + PyYAML only.
Run:  uvx piaso-mcp   (or)   python -m piaso_mcp
"""
from __future__ import annotations
import json, re, time, urllib.parse, urllib.request
from pathlib import Path
import yaml
try:                                   # MCP Python SDK 1.x
    from mcp.server.fastmcp import FastMCP as _Server
except ImportError:                    # SDK >= 2.0 (2026-07-28) renamed FastMCP -> MCPServer
    from mcp.server.mcpserver import MCPServer as _Server

DATA = Path(__file__).resolve().parent / "data"
API_BASE = "https://piaso.org/piasomarkerdb/api/v1"
REGISTRY_URL = "https://raw.githubusercontent.com/genecell/PIASO-data/master/datasets.json"
PYPI = "https://pypi.org/pypi/{name}/json"

mcp = _Server("piaso")
_cache: dict[str, tuple[float, object]] = {}


def _meta() -> dict:
    p = DATA / "meta.yaml"
    return yaml.safe_load(p.read_text()) if p.exists() else {}


def _docs() -> dict[str, str]:
    """All bundled canonical markdown, keyed by relative path."""
    out: dict[str, str] = {}
    for p in DATA.rglob("*.md"):
        out[str(p.relative_to(DATA))] = p.read_text()
    return out


def _get_json(url: str, params: dict | None = None, ttl: float = 0.0, timeout: int = 30) -> dict | list:
    q = {k: v for k, v in (params or {}).items() if v is not None}
    if q:
        url += "?" + urllib.parse.urlencode(q)
    if ttl and url in _cache and time.time() - _cache[url][0] < ttl:
        return _cache[url][1]  # type: ignore[return-value]
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "piaso-mcp"})
    with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310 (fixed hosts)
        data = json.loads(r.read().decode())
    if ttl:
        _cache[url] = (time.time(), data)
    return data


def _api_get(endpoint: str, params: dict) -> dict | list:
    return _get_json(f"{API_BASE}/{endpoint}", params)


# ------------------------------------------------------------------ doc tools
@mcp.tool()
def search_docs(query: str, max_results: int = 5) -> str:
    """Search the PIASO ecosystem knowledge pack (works with zero packages installed).

    Covers PIASO, COSG/COSGR, cytome (Python + R), LARIS, Emergene, cytorete, PIASO-data, the
    cross-component decision rules, gotchas and the piaso.org tutorial index. Returns the most
    relevant sections for `query`."""
    terms = [t.lower() for t in re.findall(r"\w+", query) if len(t) > 1]
    hits = []
    for name, text in _docs().items():
        for section in re.split(r"(?m)^(?=#{1,3}\s)", text):
            s = section.lower()
            score = sum(s.count(t) for t in terms)
            if score:
                title = section.splitlines()[0].strip() if section.strip() else name
                hits.append((score, name, title, section.strip()[:1200]))
    hits.sort(key=lambda x: -x[0])
    if not hits:
        return f"No matches for {query!r}."
    return "\n\n---\n\n".join(f"### {t}  ({n})\n{body}" for _, n, t, body in hits[:max_results])


@mcp.tool()
def get_api(function: str) -> str:
    """Data-object contract (reads/writes, defaults, what it computes) for a PIASO / COSG / cytome /
    LARIS / Emergene / cytorete function, from the component docs. The exact live signature is on
    the generated API reference at https://piaso.org/api/ (or `inspect.signature` after install)."""
    fn = function.lower().lstrip("*`")
    out = []
    for name, text in _docs().items():
        if not name.startswith("components/"):
            continue
        for section in re.split(r"(?m)^(?=#{2,3}\s)", text):
            if not section.strip():
                continue
            head = section.splitlines()[0]
            if fn in head.lower():
                out.append((0, section.strip()[:2000]))          # named in the header: rank first
            elif fn in section[:600].lower():
                out.append((1, section.strip()[:2000]))
    out.sort(key=lambda x: x[0])
    seen, uniq = set(), []
    for _, o in out:
        if o[:80] not in seen:
            seen.add(o[:80]); uniq.append(o)
    api = _meta().get("hub", {}).get("api_reference", "https://piaso.org/api/")
    if not uniq:
        return f"No API entry matching {function!r}. Live signatures: {api}"
    return "\n\n---\n\n".join(uniq) + f"\n\n(Live signature: {api})"


@mcp.tool()
def compare_implementations(function: str = "cosg") -> str:
    """COSG (Python) vs COSGR (R) divergences — params, defaults, data-object contract — and, for
    'cytome', the Python vs R cytome packages. This knowledge exists nowhere else."""
    docs = _docs()
    if function.lower().startswith("cytome"):
        text = docs.get("components/cytome.md", "")
        m = re.search(r"(?is)(## R — .*?)(\n## The Seurat|\Z)", text)
        return (m.group(1).strip()[:3500] if m else text[:3500]) or "cytome component doc not bundled."
    cosg = docs.get("components/cosg.md", "")
    if not cosg:
        return "cosg component doc not bundled."
    m = re.search(r"(?is)(## Python ↔ R divergence table.*?)(\n## Citation|\Z)", cosg)
    return (m.group(1).strip()[:3500] if m else cosg[:3500])


@mcp.tool()
def resolve_install(components: list[str], language: str = "python") -> str:
    """The exact install line(s) for a set of components in a given language (python|r), including
    extras. Answers 'I'm in R and want this chain — what do I install?' — the most error-prone thing
    in an independently-installed, cross-language ecosystem. Never suggests a matplotlib pin."""
    meta = _meta()
    comps = meta.get("components", [])
    want = {c.lower() for c in components}
    r_map = {"cosg": "cosgr", "cytome": "cytome-r"}
    lines: list[str] = []
    for c in comps:
        cid = c["id"].lower()
        keys = {cid, (c.get("pypi") or "").lower(), (c.get("name") or "").lower()}
        if language == "r":
            if "r" in c.get("language", []) and (keys & want or r_map.get(cid) in want
                                                or any(r_map.get(w) == cid for w in want)):
                lines.append(c.get("install", ""))
                for alt in c.get("install_alternatives", []):
                    lines.append(f"# alternative: {alt}")
        elif "python" in c.get("language", []) and keys & want:
            lines.append(c.get("install", ""))
            for k, v in (c.get("extras") or {}).items():
                lines.append(f"# extra [{k}]: {v}")
    if language == "r" and (want & {"piaso", "laris", "emergene", "cytorete", "scalar", "gdr", "infog"}):
        lines.append("# INFOG / GDR / PIASOscore / SCALAR / LARIS / cytorete are Python-only: write a .cytome from R "
                     "(write_cytome(obj, path)), run the Python call on it, read results back with read_cytome().")
    if not lines:
        return f"No {language} install found for {components}. Known: {[c['id'] for c in comps]}"
    seen: list[str] = []
    for ln in lines:
        if ln and ln not in seen:
            seen.append(ln)
    return "\n".join(seen)


@mcp.tool()
def list_tutorials(topic: str | None = None, component: str | None = None) -> str:
    """The piaso.org tutorial index (executed, human-reviewed): title, URL, what it covers, which
    components it uses. Filter by topic (scRNA-seq, methods, marker-genes, annotation, gene-sets,
    spatial, grn, cell-cell-interaction, plotting-data) and/or component (piaso, cosg, cytome,
    laris, cytorete, emergene). Route the user to the tutorial for their platform before writing code."""
    meta = _meta()
    base = meta.get("hub", {}).get("tutorials", "https://piaso.org/tutorials/").rstrip("/") + "/"
    rows = []
    for t in meta.get("tutorials", []):
        if topic and t.get("topic", "").lower() != topic.lower():
            continue
        if component and component.lower() not in [c.lower() for c in t.get("components", [])]:
            continue
        rows.append(f"- [{t['topic']}] {t['title']} — {base}{t['slug']}/\n    {t['covers']} ({', '.join(t.get('components', []))})")
    if not rows:
        return f"No tutorials for topic={topic!r} component={component!r}. Index: {base}"
    return "\n".join(rows)


@mcp.tool()
def version_matrix() -> str:
    """The component versions this knowledge pack was tested against (from meta.yaml), with install
    lines, roles and citations — so an agent can state what the docs assume."""
    meta = _meta()
    hub = meta.get("hub", {})
    out = [f"{hub.get('name')} {hub.get('version')} — tested on {hub.get('tested_on')}"]
    for c in meta.get("components", []):
        out.append(f"- {c['name']} ({c['id']}): {c.get('version_last_tested')} | {'/'.join(c.get('language', []))} | "
                   f"install `{c.get('install')}` | cite {c.get('citation_id')}")
    return "\n".join(out)


@mcp.tool()
def check_versions() -> str:
    """Current PyPI release of each Python package vs the version this snapshot was tested against.
    Pure PyPI JSON (network); use it to disclose drift ('docs tested on 1.2.3, you have 1.3.0')."""
    meta = _meta()
    rows = []
    for c in meta.get("components", []):
        name = c.get("pypi")
        if not name:
            continue
        try:
            latest = _get_json(PYPI.format(name=name), ttl=3600)["info"]["version"]  # type: ignore[index]
        except Exception as e:  # noqa: BLE001
            latest = f"unavailable ({type(e).__name__})"
        tested = c.get("version_last_tested")
        flag = "" if latest == tested else "  <- differs from tested version"
        rows.append(f"- {name}: PyPI {latest} | tested {tested}{flag}")
    return "\n".join(rows) if rows else "No PyPI components in meta.yaml."


# ------------------------------------------------------------- registry tools
def _registry() -> dict:
    return _get_json(REGISTRY_URL, ttl=86400)  # type: ignore[return-value]


@mcp.tool()
def list_datasets() -> str:
    """PIASO-data datasets (live registry from genecell/PIASO-data, cached 24 h): id, format, size,
    species, cells, counts layer, how to load. Falls back to the bundled data.md offline."""
    try:
        reg = _registry()
        items = reg.get("datasets", [])
        if isinstance(items, dict):
            items = [dict(id=k, **v) for k, v in items.items()]
        lines = [f"PIASO-data registry v{reg.get('version')} — {reg.get('zenodo_record')} (concept DOI {reg.get('zenodo_doi')})",
                 "load: piaso.data.load_dataset(id)  |  piaso.data.load_dataset(id, return_type='cytome')  |  piaso.data.fetch_dataset(id)"]
        for x in items:
            mb = (x.get("size_bytes") or 0) / 1e6
            lines.append(f"- {x.get('id')}: {x.get('format')} | {mb:,.1f} MB | {x.get('species')} | cells={x.get('cells')} | "
                         f"counts_layer={x.get('counts_layer')} | {x.get('title')}")
        return "\n".join(lines)
    except Exception as e:  # noqa: BLE001
        data_md = _docs().get("data.md", "")
        return f"(live registry unavailable: {e}; bundled snapshot follows)\n\n" + data_md[:3500]


@mcp.tool()
def get_dataset(name: str) -> str:
    """One PIASO-data registry entry (live): title, url, md5, size, cells, features, counts layer,
    reference, tutorials that use it."""
    try:
        reg = _registry()
        items = reg.get("datasets", [])
        if isinstance(items, dict):
            items = [dict(id=k, **v) for k, v in items.items()]
        for x in items:
            if x.get("id") == name:
                return json.dumps(x, indent=1)
        return f"No dataset {name!r}. Known: {[x.get('id') for x in items]}"
    except Exception as e:  # noqa: BLE001
        return f"live registry unavailable: {e}"


# --------------------------------------------------------------- marker tools
@mcp.tool()
def query_marker_db(gene: str | None = None, cell_type: str | None = None,
                    study: str | None = None, species: str | None = None,
                    tissue: str | None = None, limit: int = 20) -> str:
    """Query the live PIASOmarkerDB for cell-type marker genes (proxies the piaso.org REST API).

    Works with zero packages installed, and is the programmatic path for R / non-Python users.
    Returns rows with: cell_type, condition, gene, species, specificity_score (a COSG score),
    study_publication, tissue. Cell-type names must match the study's vocabulary exactly.
    Requires network access to piaso.org."""
    try:
        rows = _api_get("markers", dict(gene=gene, cell_type=cell_type, study=study,
                                        species=species, tissue=tissue, limit=limit))
    except Exception as e:  # noqa: BLE001
        return f"PIASOmarkerDB request failed: {e}"
    return json.dumps(rows, indent=1)[:4000]


@mcp.tool()
def get_markers(gene: str | None = None, cell_type: str | None = None,
                species: str | None = None, limit: int = 20) -> str:
    """Alias of query_marker_db (mirrors piaso.tl.getMarkers)."""
    return query_marker_db(gene=gene, cell_type=cell_type, species=species, limit=limit)


@mcp.tool()
def list_studies() -> str:
    """List the studies available in PIASOmarkerDB (proxies the live API; 36 at last check)."""
    try:
        return json.dumps(_api_get("studies", {}), indent=1)[:4000]
    except Exception as e:  # noqa: BLE001
        return f"PIASOmarkerDB request failed: {e}"


def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()
