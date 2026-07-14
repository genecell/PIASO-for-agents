"""piaso-mcp — a local stdio MCP server for the PIASO ecosystem.

Serves KNOWLEDGE and PUBLIC DATA only — never compute, never user data.
- Doc tools (search_docs, get_api, compare_implementations, resolve_install,
  list_datasets) read a bundled snapshot of canonical/ (populated by build.py).
- Marker tools (query_marker_db, get_markers, list_studies) PROXY the live
  PIASOmarkerDB REST API at https://piaso.org/piasomarkerdb/api/v1/ (decision:
  proxy the existing API rather than bundle the data). This gives R / non-Python
  callers their first programmatic path to PIASOmarkerDB.

Needs NO PIASO packages installed. Stdlib + the `mcp` SDK + PyYAML only.
Run:  uvx piaso-mcp   (or)   python -m piaso_mcp
"""
from __future__ import annotations
import json, re, urllib.parse, urllib.request
from pathlib import Path
import yaml
from mcp.server.fastmcp import FastMCP

DATA = Path(__file__).resolve().parent / "data"
API_BASE = "https://piaso.org/piasomarkerdb/api/v1"

mcp = FastMCP("piaso")


def _meta() -> dict:
    p = DATA / "meta.yaml"
    return yaml.safe_load(p.read_text()) if p.exists() else {}


def _docs() -> dict[str, str]:
    """All bundled canonical markdown, keyed by relative path."""
    out: dict[str, str] = {}
    for p in DATA.rglob("*.md"):
        out[str(p.relative_to(DATA))] = p.read_text()
    return out


def _api_get(endpoint: str, params: dict) -> dict:
    q = {k: v for k, v in params.items() if v is not None}
    url = f"{API_BASE}/{endpoint}"
    if q:
        url += "?" + urllib.parse.urlencode(q)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:  # noqa: S310 (fixed host)
        return json.loads(r.read().decode())


# ------------------------------------------------------------------ doc tools
@mcp.tool()
def search_docs(query: str, max_results: int = 5) -> str:
    """Search the PIASO ecosystem knowledge pack (works with zero packages installed).

    Returns the most relevant sections from the canonical docs for `query`."""
    terms = [t.lower() for t in re.findall(r"\w+", query) if len(t) > 1]
    hits = []
    for name, text in _docs().items():
        # split into sections by markdown headers
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
    """Data-object contract (reads/writes, what it computes) for a PIASO/COSG/LARIS/Emergene
    function, from the component docs. For the exact signature, install the package and use
    `inspect.signature` — this returns the documented pre/postconditions."""
    fn = function.lower().lstrip("*`")
    out = []
    for name, text in _docs().items():
        if not name.startswith("components/"):
            continue
        for section in re.split(r"(?m)^(?=#{2,3}\s)", text):
            if not section.strip():
                continue
            head = section.splitlines()[0]
            # match the header (function name in the title) or an early body mention
            if fn in head.lower() or fn in section[:600].lower():
                out.append(section.strip()[:2000])
    # de-dup while preserving order
    seen, uniq = set(), []
    for o in out:
        if o[:80] not in seen:
            seen.add(o[:80]); uniq.append(o)
    return "\n\n---\n\n".join(uniq) if uniq else f"No API entry matching {function!r}."


@mcp.tool()
def compare_implementations(function: str = "cosg") -> str:
    """COSG (Python) vs COSGR (R) divergences — params, defaults, data-object contract.

    This knowledge exists nowhere else: the two implementations have different defaults."""
    cosg = _docs().get("components/cosg.md", "")
    if not cosg:
        return "cosg component doc not bundled."
    # return the divergence table + surrounding context
    m = re.search(r"(?is)(divergen.*?)(\n#{1,2}\s|\Z)", cosg)
    return (m.group(1).strip()[:3000] if m else cosg[:3000])


@mcp.tool()
def resolve_install(components: list[str], language: str = "python") -> str:
    """The exact install line for a set of components in a given language (python|r).

    Answers 'I'm in R and want this chain — what do I install?' — the most error-prone
    thing in an independently-installed, cross-language ecosystem."""
    meta = _meta()
    by_id = {c["id"]: c for c in meta.get("components", [])}
    want = [c.lower() for c in components]
    lines: list[str] = []
    for cid, c in by_id.items():
        if cid in want or c.get("pypi", "") in want:
            langs = c.get("language", [])
            if language == "r" and "r" in langs:
                lines.append(c.get("install", ""))
            elif language != "r" and "python" in langs:
                lines.append(c.get("install", ""))
    # handle cosg/cosgr language routing
    if language == "r" and any(w in ("cosg", "cosgr") for w in want):
        lines.append('remotes::install_github("genecell/COSGR")')
    if not lines:
        return f"No {language} install found for {components}. Known: {list(by_id)}"
    seen = []
    for ln in lines:
        if ln and ln not in seen:
            seen.append(ln)
    return "\n".join(seen)


@mcp.tool()
def list_datasets() -> str:
    """PIASO-data fixtures (URLs + loading notes) so an agent can write a runnable example."""
    data_md = _docs().get("data.md", "")
    return data_md[:3000] if data_md else "data.md not bundled."


# --------------------------------------------------------------- marker tools
@mcp.tool()
def query_marker_db(gene: str | None = None, cell_type: str | None = None,
                    study: str | None = None, species: str | None = None,
                    tissue: str | None = None, limit: int = 20) -> str:
    """Query the live PIASOmarkerDB for cell-type marker genes (proxies piaso.org API).

    Works with zero packages installed, and is the first programmatic path for R / non-Python
    users. Returns rows with: cell_type, gene, species, specificity_score, study_publication,
    tissue. Requires network access to piaso.org."""
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
    """List the studies available in PIASOmarkerDB (proxies the live API)."""
    try:
        return json.dumps(_api_get("studies", {}), indent=1)[:4000]
    except Exception as e:  # noqa: BLE001
        return f"PIASOmarkerDB request failed: {e}"


def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()
