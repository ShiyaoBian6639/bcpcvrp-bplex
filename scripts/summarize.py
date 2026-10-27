#!/usr/bin/env python3
"""Summarize cvrpbcp results.jsonl: table to stdout, proven-only usable.csv.

    python -I scripts/summarize.py results/A.jsonl [more.jsonl ...]

A run is "usable" when status is optimal, the fleet cap was not binding and,
if a published value is known, the proven cost equals it.
"""
import csv
import json
import sys


def main(paths):
    rows = {}
    for p in paths:
        with open(p, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                r = json.loads(line)
                rows[r["instance"]] = r  # last run wins
    rows = sorted(rows.values(), key=lambda r: (r["n"], r["instance"]))
    print(f"{'instance':14} {'n':>4} {'status':10} {'lb':>8} {'ub':>8} {'pub':>8} {'nodes':>6} {'time':>8} {'root':>8} verdict")
    usable = []
    nopt = nmatch = 0
    for r in rows:
        pub = r.get("published", -1)
        verdict = ""
        if r["status"] == "optimal":
            nopt += 1
            if pub >= 0:
                verdict = "match" if r["ub"] == pub else "DIFFERS"
                nmatch += r["ub"] == pub
            if r.get("fleet_cap_binding"):
                verdict += " CAP-BINDING"
            if (pub < 0 or r["ub"] == pub) and not r.get("fleet_cap_binding"):
                usable.append(r)
        elif pub >= 0 and r["lb"] > pub:
            verdict = "LB>PUB BUG"
        print(f"{r['instance']:14} {r['n']:4d} {r['status']:10} {r['lb']:8d} {r['ub']:8d} {pub:8d} {r['nodes']:6d} "
              f"{r['time']:8.1f} {r['time_root']:8.1f} {verdict}")
    print(f"\n{len(rows)} runs, {nopt} optimal, {nmatch} match published, {len(usable)} usable")
    out = paths[0].rsplit(".", 1)[0] + ".usable.csv"
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["instance", "n", "cost", "scale", "nodes", "time", "root_lb"])
        for r in usable:
            w.writerow([r["instance"], r["n"], r["ub"], r.get("scale", 1), r["nodes"], f"{r['time']:.2f}", f"{r['root_lb']:.4f}"])
    print(f"wrote {out}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    main(sys.argv[1:])
