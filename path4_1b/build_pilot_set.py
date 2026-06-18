#!/usr/bin/env python3
"""Blended pilot set: consensus + per-axis picks + model/tier diversity.

The per-model-z and ESM-only shortlists overlap only 4/20 — the scoring choice is
consequential and predictor-vs-function is unvalidated at this divergence. So the
pilot is designed to LEARN which axis tracks function: it samples both philosophies
plus their consensus, with deliberate model-type and tier (bold/safe) spread.

Selection (reproducible):
  * consensus   : designs in BOTH shortlists (safest; win either way)
  * esm-axis    : from ESM-only-unique, top by ESM, prioritising Soluble/Protein
                  (the variants the biased MPNN axis hid)
  * z-axis      : from per-model-z-unique, top by combined z, both tiers
Writes pilot_set.json ({"shortlist":[...]} so build_af3.py can read it) + audits all.

Usage: python build_pilot_set.py [--n-esm 5 --n-z 5]
"""
import argparse, csv, json, os
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))


def diverse_pick(cands, n, key, min_nonligand=2):
    """top-n by key, but greedily ensure >= min_nonligand non-Ligand if available."""
    cands = sorted(cands, key=lambda r: -key(r))
    chosen, nonlig = [], 0
    # first pass: take non-Ligand to meet the floor
    for r in cands:
        if len(chosen) >= n: break
        if r["model"] != "Ligand" and nonlig < min_nonligand:
            chosen.append(r); nonlig += 1
    # fill remaining by score
    for r in cands:
        if len(chosen) >= n: break
        if r not in chosen:
            chosen.append(r)
    return chosen


def main(a):
    rows = {r["id"]: r for r in csv.DictReader(open(os.path.join(HERE, "pareto_ranked.csv")))}
    for r in rows.values():
        r["esm"] = float(r["esm"]); r["z"] = float(r["z"]); r["pct_identity"] = float(r["pct_identity"])
    pmz = json.load(open(os.path.join(HERE, "synthesis_shortlist.json")))["shortlist"]
    esm = json.load(open(os.path.join(HERE, "esm_only_shortlist.json")))["shortlist"]
    pmz_ids = {d["id"] for d in pmz}; esm_ids = {d["id"] for d in esm}
    consensus = sorted(pmz_ids & esm_ids)
    esm_unique = [rows[i] for i in (esm_ids - pmz_ids)]
    z_unique = [rows[i] for i in (pmz_ids - esm_ids)]

    picks = {}
    for i in consensus:
        picks[i] = "consensus"
    for r in diverse_pick(esm_unique, a.n_esm, lambda r: r["esm"]):
        picks[r["id"]] = "esm-axis"
    for r in diverse_pick(z_unique, a.n_z, lambda r: r["z"]):
        picks.setdefault(r["id"], "z-axis")
    # bold-novelty bucket: ensure the high-novelty T_aggressive flavor (~57-61% id, riskiest
    # + the patent story) is sampled. Lowest-identity frontier designs not already chosen.
    bold_pool = [r for r in (list(esm_unique) + list(z_unique)) + [rows[i] for i in consensus]
                 if r["tier"] == "T_aggressive" and r["id"] not in picks]
    for r in sorted(bold_pool, key=lambda r: r["pct_identity"])[:a.n_bold]:
        picks[r["id"]] = "bold-novelty"

    meta = {r["id"]: r for r in json.load(open(os.path.join(HERE, "unique_meta.json")))}
    r2i = {int(k): v for k, v in json.load(open(os.path.join(HERE, "wt_resnum_to_index.json"))).items()}
    CAT = {10:'D',762:'E',839:'D',840:'H',863:'N',983:'H',986:'D'}
    SPRY = {61:'R',1111:'R',1135:'L',1136:'W',1218:'K',1219:'Q',1317:'R',1322:'R',1333:'P',1335:'Q',1337:'R'}
    out, bad = [], 0
    for sid, prov in picks.items():
        seq = meta[sid]["seq"]; r = rows[sid]
        ok = (len(seq) == 1341
              and all(seq[r2i[p]] == aa for p, aa in CAT.items())
              and all(seq[r2i[p]] == aa for p, aa in SPRY.items())
              and [p for p, i in r2i.items() if seq[i] == 'C'] == [574])
        if not ok:
            bad += 1; print("AUDIT FAIL", sid)
        out.append(dict(id=sid, provenance=prov, tier=r["tier"], model=r["model"],
                        esm=r["esm"], pct_identity=r["pct_identity"], seq=seq))
    out.sort(key=lambda d: (d["provenance"], -d["esm"]))
    print(f"pilot set: {len(out)} designs")
    print("  by provenance:", dict(Counter(d['provenance'] for d in out)))
    print("  by model:", dict(Counter(d['model'] for d in out)))
    print("  by tier:", dict(Counter(d['tier'] for d in out)))
    print(f"  identity range: {min(d['pct_identity'] for d in out):.1f}-{max(d['pct_identity'] for d in out):.1f}%")
    print(f"  sequence audit: {'ALL PASS' if bad == 0 else str(bad)+' FAIL'} (len 1341, 7/7 catalytic, 11/11 SpRY, omit-C)")
    json.dump({"description": "blended pilot: consensus + esm-axis + z-axis (learns predictor-vs-function)",
               "shortlist": out}, open(os.path.join(HERE, "pilot_set.json"), "w"), indent=2)
    print("wrote pilot_set.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-esm", type=int, default=5, dest="n_esm")
    ap.add_argument("--n-z", type=int, default=5, dest="n_z")
    ap.add_argument("--n-bold", type=int, default=4, dest="n_bold")
    main(ap.parse_args())
