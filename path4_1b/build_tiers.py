#!/usr/bin/env python3
"""Per-domain graded tiers for the SpRY-Cas9 Path-4.1b sweep.

Implements the chosen "per-domain identity" strategy by applying a *different*
conservation floor to each domain (PROTECT vs MODERATE vs REDESIGN from
domain_map.json), then layering the always-fixed functional sets:

    fixed(tier) =  per-domain-conservation(stringency)
                 U SpRY(11)  U catalytic(7)  U contacts(35)  U DCA(42)

We emit a graded SWEEP of 3 tiers (T_aggressive / T_balanced / T_safe) by
shifting every domain's floor up/down one notch. T_balanced is the recommended
operating point (mirrors how RT-4.1b's +50 tier was the sweet spot, but here the
"+50" is applied only to the PROTECT domains, not globally).

Inputs (all already author-numbered, no MSA re-run):
  domain_map.json                      (this dir; from domain_map.py)
  ../outputs/constraints/8srs_constraints.json   (conservation_tiers, dca_top, spry)
  ../fixed_residues/tier1_fixed.txt    (to recover the 35 nucleic-acid/active-site contacts)

Outputs:
  tiers/<tier>_fixed.txt   LigandMPNN fixed string ("A10 A56 ...")
  tiers_summary.json       per-domain + overall fixed counts / fixed-fraction
Fixed-fraction is the pre-generation proxy for identity-to-WT (fixed positions
stay WT); ACTUAL %id is measured after MPNN, as in the RT campaign.
"""
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

dm = json.load(open(os.path.join(HERE, "domain_map.json")))
con = json.load(open(os.path.join(REPO, "outputs/constraints/8srs_constraints.json")))

SPRY = set(p["pos"] for p in dm["spry_positions"])
CATALYTIC = set(c["pos"] for c in dm["catalytic"])
DCA = set(con["dca_top"])
# conservation sets (author resnums) keyed by integer floor percent
CONS = {90: set(con["conservation_tiers"]["ge_90"]),
        70: set(con["conservation_tiers"]["ge_70"]),
        50: set(con["conservation_tiers"]["ge_50"])}

# the 35 nucleic-acid / active-site contacts = tier1 minus SpRY minus catalytic
tier1 = {int(tok[1:]) for tok in open(os.path.join(REPO, "fixed_residues/tier1_fixed.txt")).read().split()}
CONTACTS = tier1 - SPRY - CATALYTIC

# domain -> (start, end, policy)
DOMS = [(d["name"], d["start"], d["end"], d["policy"]) for d in dm["domains"]]

# per-tier conservation floor by policy. Lower floor = MORE fixed = safer/less novel.
# stringency shifts all three together. floor 100 == "fix nothing by conservation".
TIER_FLOORS = {
    # tier name        PROTECT  MODERATE  REDESIGN
    "T_aggressive": dict(PROTECT=70, MODERATE=90, REDESIGN=100),  # max novelty
    "T_balanced":   dict(PROTECT=50, MODERATE=70, REDESIGN=90),   # <-- operating point
    "T_safe":       dict(PROTECT=50, MODERATE=50, REDESIGN=70),   # most likely functional
}

ALWAYS = SPRY | CATALYTIC | CONTACTS | DCA  # functional locks, every tier


def cons_at(floor):
    """conserved set at a given floor; floor==100 -> empty (fix nothing by conservation)."""
    return set() if floor >= 100 else CONS[floor]


def in_domain(pos, lo, hi):
    return lo <= pos <= hi


def build_tier(floors):
    fixed = set(ALWAYS)
    per_dom = {}
    for name, lo, hi, pol in DOMS:
        floor = floors[pol]
        dom_cons = {p for p in cons_at(floor) if in_domain(p, lo, hi)}
        fixed |= dom_cons
        dom_len = hi - lo + 1
        dom_fixed = {p for p in fixed if in_domain(p, lo, hi)}
        per_dom[name] = dict(policy=pol, floor=floor, length=dom_len,
                             fixed=len(dom_fixed),
                             fixed_frac=round(len(dom_fixed) / dom_len, 3))
    return fixed, per_dom


def main():
    os.makedirs(os.path.join(HERE, "tiers"), exist_ok=True)
    L = dm["full_length_ref"]
    summary = {"target": dm["target"], "always_fixed": dict(
        spry=len(SPRY), catalytic=len(CATALYTIC), contacts=len(CONTACTS),
        dca=len(DCA), union=len(ALWAYS)), "tiers": {}}
    print(f"always-fixed functional locks: SpRY {len(SPRY)} + catalytic {len(CATALYTIC)} "
          f"+ contacts {len(CONTACTS)} + DCA {len(DCA)} = {len(ALWAYS)} union\n")
    for tname, floors in TIER_FLOORS.items():
        fixed, per_dom = build_tier(floors)
        frac = len(fixed) / L
        # emit LigandMPNN fixed string
        s = " ".join(f"A{p}" for p in sorted(fixed))
        open(os.path.join(HERE, "tiers", f"{tname}_fixed.txt"), "w").write(s + "\n")
        summary["tiers"][tname] = dict(floors=floors, n_fixed=len(fixed),
                                       fixed_frac=round(frac, 3), per_domain=per_dom)
        print(f"=== {tname}  ({len(fixed)} fixed, {100*frac:.1f}% of {L}) ===")
        print(f"  floors: PROTECT>={floors['PROTECT']} MODERATE>={floors['MODERATE']} "
              f"REDESIGN>={floors['REDESIGN']}")
        for name, lo, hi, pol in DOMS:
            d = per_dom[name]
            print(f"    {name:12} {pol:9} {d['fixed']:>4}/{d['length']:<4} "
                  f"({100*d['fixed_frac']:>4.0f}% fixed)")
        print()
    json.dump(summary, open(os.path.join(HERE, "tiers_summary.json"), "w"), indent=2)
    # sanity: nesting of functional locks; every SpRY+catalytic present in every tier
    for tname in TIER_FLOORS:
        fx = {int(t[1:]) for t in open(os.path.join(HERE, "tiers", f"{tname}_fixed.txt")).read().split()}
        assert SPRY <= fx and CATALYTIC <= fx, f"{tname} dropped a functional lock!"
    print("AUDIT: all tiers retain 11 SpRY + 7 catalytic. wrote tiers/ + tiers_summary.json")


if __name__ == "__main__":
    main()
