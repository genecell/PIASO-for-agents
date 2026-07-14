"""Score blind-router verdicts against tests/activation_prompts.yaml.

Usage: python tests/score_activation.py judgeA.json judgeB.json judgeC.json
Each judge file is a JSON array of {id, activate, component, language, reason}.
Prints a per-prompt table (expected vs majority vote), accuracy, and any mismatch.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import yaml

HERE = Path(__file__).resolve().parent
SPEC = {p["id"]: p for p in yaml.safe_load((HERE / "activation_prompts.yaml").read_text())["prompts"]}


def norm(x):
    return (x or "").strip().lower() if isinstance(x, str) else x


def load(paths):
    judges = []
    for p in paths:
        arr = json.loads(Path(p).read_text())
        judges.append({d["id"]: d for d in arr})
    return judges


def majority(vals):
    t = sum(1 for v in vals if v)
    return t > len(vals) / 2, t, len(vals)


def main(paths):
    judges = load(paths)
    fire_ok = route_ok = lang_ok = 0
    fire_total = route_total = lang_total = 0
    rows = []
    for pid, spec in SPEC.items():
        verdicts = [j.get(pid, {}) for j in judges]
        acts = [bool(v.get("activate")) for v in verdicts]
        maj_fire, nyes, n = majority(acts)
        exp_fire = spec["fire"]
        fire_total += 1
        f_ok = maj_fire == exp_fire
        fire_ok += f_ok
        # component/language only checked when it should fire
        r_ok = l_ok = None
        if exp_fire:
            exp_comp = [c.lower() for c in spec.get("component", [])]
            comps = [norm(v.get("component")) for v in verdicts if v.get("activate")]
            # majority component vote among activating judges
            r_ok = bool(comps) and any(c in exp_comp for c in _mode(comps))
            route_total += 1; route_ok += bool(r_ok)
            if "language" in spec:
                langs = [norm(v.get("language")) for v in verdicts if v.get("activate")]
                l_ok = bool(langs) and _mode(langs)[0] == norm(spec["language"])
                lang_total += 1; lang_ok += bool(l_ok)
        rows.append((pid, exp_fire, f"{nyes}/{n}", maj_fire, f_ok, r_ok, l_ok,
                     _mode([norm(v.get("component")) for v in verdicts if v.get("activate")])))
    # report
    print(f"{'prompt':22} {'exp':4} {'votes':6} {'fire?':5} {'route':6} {'lang':5}  component(maj)")
    print("-" * 78)
    for pid, ef, votes, mf, fok, rok, lok, comp in rows:
        mark = lambda b: "·" if b is None else ("✓" if b else "✗")
        print(f"{pid:22} {str(ef):4} {votes:6} {mark(fok):5} {mark(rok):6} {mark(lok):5}  {comp}")
    print("-" * 78)
    print(f"Activation (fire/no-fire): {fire_ok}/{fire_total}")
    print(f"Component routing (when firing): {route_ok}/{route_total}")
    print(f"Language routing (probes): {lang_ok}/{lang_total}")
    ok = fire_ok == fire_total and route_ok == route_total and lang_ok == lang_total
    print("RESULT:", "ALL PASS" if ok else "MISMATCHES ABOVE")
    return 0 if ok else 1


def _mode(xs):
    xs = [x for x in xs if x]
    if not xs:
        return [None]
    from collections import Counter
    c = Counter(xs)
    top = c.most_common(1)[0][1]
    return [k for k, v in c.items() if v == top]


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
