#!/usr/bin/env python3
"""ESM-only cross-check shortlist (single neutral judge) vs the per-model-z shortlist.

The per-model-z fix patches the cross-model MPNN-score incomparability statistically.
This is the clean alternative: rank ALL designs by the model-AGNOSTIC ESM naturalness
axis alone (one consistent judge for every design), apply the <70% identity patent
gate, take the top N. Then compare overlap with synthesis_shortlist.json — high
overlap = the 20-design pick is robust to the scoring choice; the intersection is the
safest set to carry forward.

Reads pareto_ranked.csv (already has esm + pct_identity + model for all 1530) so no
recompute. Writes esm_only_shortlist.json + overlap report. New (non-overlapping)
designs are sequence-audited.

Usage: python esm_only_shortlist.py [--top 20 --id-ceiling 70]
"""
import argparse, csv, json, os
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))


def main(a):
    rows = list(csv.DictReader(open(os.path.join(HERE, "pareto_ranked.csv"))))
    for r in rows:
        r["esm"] = float(r["esm"]); r["pct_identity"] = float(r["pct_identity"])
    gated = [r for r in rows if r["pct_identity"] < a.id_ceiling]
    gated.sort(key=lambda r: -r["esm"])              # most natural first
    top = gated[:a.top]

    pmz = json.load(open(os.path.join(HERE, "synthesis_shortlist.json")))["shortlist"]
    pmz_ids = {d["id"] for d in pmz}
    esm_ids = {r["id"] for r in top}
    overlap = esm_ids & pmz_ids

    print(f"ESM-only top-{a.top} (<{a.id_ceiling}% id) by MODEL: {dict(Counter(r['model'] for r in top))}")
    print(f"per-model-z shortlist ({len(pmz)}) by MODEL: {dict(Counter(d['model'] for d in pmz))}")
    print(f"\nOVERLAP: {len(overlap)}/{a.top} designs shared between the two shortlists")
    print(f"  shared ids: {sorted(overlap)}")
    print(f"  ESM-only ONLY: {sorted(esm_ids - pmz_ids)}")
    print(f"  per-model-z ONLY: {sorted(pmz_ids - esm_ids)}")

    # write + audit the ESM-only shortlist
    meta = {r["id"]: r for r in json.load(open(os.path.join(HERE, "unique_meta.json")))}
    r2i = {int(k): v for k, v in json.load(open(os.path.join(HERE, "wt_resnum_to_index.json"))).items()}
    CAT = {10:'D',762:'E',839:'D',840:'H',863:'N',983:'H',986:'D'}
    SPRY = {61:'R',1111:'R',1135:'L',1136:'W',1218:'K',1219:'Q',1317:'R',1322:'R',1333:'P',1335:'Q',1337:'R'}
    bad = 0; out = []
    for r in top:
        seq = meta[r["id"]]["seq"]
        ok = (len(seq) == 1341
              and all(seq[r2i[p]] == aa for p, aa in CAT.items())
              and all(seq[r2i[p]] == aa for p, aa in SPRY.items())
              and [p for p, i in r2i.items() if seq[i] == 'C'] == [574])
        if not ok:
            bad += 1; print("  AUDIT FAIL", r["id"])
        out.append(dict(id=r["id"], tier=r["tier"], model=r["model"],
                        esm=r["esm"], pct_identity=r["pct_identity"],
                        in_per_model_z=(r["id"] in pmz_ids), seq=seq))
    print(f"\nsequence audit of ESM-only top-{a.top}: {'ALL PASS' if bad == 0 else str(bad)+' FAIL'} "
          "(len 1341, 7/7 catalytic, 11/11 SpRY, omit-C)")
    json.dump({"axis": "ESM naturalness only (model-agnostic single judge)",
               "id_ceiling": a.id_ceiling, "overlap_with_per_model_z": len(overlap),
               "intersection_ids": sorted(overlap), "shortlist": out},
              open(os.path.join(HERE, "esm_only_shortlist.json"), "w"), indent=2)
    print("wrote esm_only_shortlist.json")
    print(f"\nROBUSTNESS: {len(overlap)}/{a.top} agreement -> intersection ({len(overlap)} designs) is the "
          "safest carry-forward set; union covers both scoring philosophies.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--id-ceiling", type=float, default=70.0, dest="id_ceiling")
    main(ap.parse_args())
