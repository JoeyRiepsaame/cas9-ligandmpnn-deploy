#!/usr/bin/env python3
"""Pre-design blind-spot audit for the SpRY-Cas9 Path-4.1b tiers.

Cheap insurance BEFORE spending GPU: predict which native residues the
generation constraints will forcibly mutate, and whether that creates artificial
co-variation (the RT campaign's omit-C / E-penalty artifact failure mode).

Path-4.1b constraints here:
  --omit_AA C      -> Cat1 (omit) artifacts possible
  NO --bias_AA     -> Cat2 (bias) artifacts ELIMINATED by design (E-penalty dropped)

For each tier we:
  1. classify every native cysteine as FIXED (protected) or FREE (forced away by omit-C),
  2. run the universal pre_design_scan (forced-departure victims + predicted correlations),
  3. cross-check each forced cysteine against Cas9-family conservation (is evolution
     telling us to keep it?) and its domain,
  4. confirm bias artifacts == 0 (validates the E-penalty drop).
Recommendation: lock any FREE cysteine that is conserved (>=50% family) before generation.
"""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from mpnn_blind_spot_detector import pre_design_scan

OMIT = "C"
BIAS = None   # Path-4.1b: E-penalty dropped

wt = "".join(l.strip() for l in open(os.path.join(HERE, "spry_cas9_wt.fasta")) if not l.startswith(">"))
r2i = {int(k): v for k, v in json.load(open(os.path.join(HERE, "wt_resnum_to_index.json"))).items()}
i2r = {v: k for k, v in r2i.items()}
dm = json.load(open(os.path.join(HERE, "domain_map.json")))
con = json.load(open(os.path.join(REPO, "outputs/constraints/8srs_constraints.json")))

# conservation status per author resnum
CONS = {90: set(con["conservation_tiers"]["ge_90"]),
        70: set(con["conservation_tiers"]["ge_70"]),
        50: set(con["conservation_tiers"]["ge_50"])}
def cons_label(resnum):
    for f in (90, 70, 50):
        if resnum in CONS[f]:
            return f">={f}%"
    return "<50%"
def domain_of(resnum):
    for d in dm["domains"]:
        if d["start"] <= resnum <= d["end"]:
            return d["name"], d["policy"]
    return "?", "?"

native_cys = sorted(i2r[i] for i, a in enumerate(wt) if a == "C")
print(f"WT SpRY-Cas9: {len(wt)} aa, native cysteines = {native_cys} "
      f"(only {len(native_cys)} -> omit-C artifact surface is structurally small)\n")
for c in native_cys:
    dom, pol = domain_of(c)
    print(f"  C{c}: domain {dom} ({pol}), family conservation {cons_label(c)}")

TIERS = ["T_aggressive", "T_balanced", "T_safe"]
summary = {"omit": OMIT, "bias": BIAS, "native_cys": native_cys, "tiers": {}}
concerns = []
for t in TIERS:
    fixed_resnums = {int(tok[1:]) for tok in open(os.path.join(HERE, "tiers", f"{t}_fixed.txt")).read().split()}
    # pre_design_scan wants 1-indexed WT-sequence positions
    fixed_1idx = sorted(r2i[r] + 1 for r in fixed_resnums if r in r2i)
    scan = pre_design_scan(wt, fixed_1idx, omit_aas=OMIT, bias_aas=BIAS)
    # map forced-departure warnings (seq 1-idx) back to author resnums
    forced = []
    for w in scan["warnings"]:
        rn = i2r[w["position_1idx"] - 1]
        forced.append(rn)
    free_cys = [c for c in native_cys if c not in fixed_resnums]
    prot_cys = [c for c in native_cys if c in fixed_resnums]
    print(f"\n=== {t} ===")
    print(f"  cysteines: protected(fixed)={prot_cys}  free(forced by omit-C)={free_cys}")
    print(f"  forced-departure victims (should equal free cys): {sorted(forced)}")
    print(f"  predicted artificial correlations: {scan['total_predicted_correlations']}  "
          f"| bias artifacts: {sum(1 for w in scan['warnings'] if w['impact']=='PENALIZED')} (expect 0)")
    for c in free_cys:
        lab = cons_label(c); dom, pol = domain_of(c)
        risky = lab != "<50%"
        flag = "  <-- CONSIDER LOCKING (conserved+free)" if risky else ""
        print(f"    free C{c}: {dom}/{pol}, conservation {lab}{flag}")
        if risky:
            concerns.append((t, c, lab))
    # sanity: forced victims are exactly the free cysteines (omit C only)
    assert set(forced) == set(free_cys), f"{t}: forced set {sorted(forced)} != free cys {free_cys}"
    assert sum(1 for w in scan["warnings"] if w["impact"] == "PENALIZED") == 0, "unexpected bias artifact"
    summary["tiers"][t] = dict(protected_cys=prot_cys, free_cys=free_cys,
                               forced=sorted(forced),
                               predicted_correlations=scan["total_predicted_correlations"])

json.dump(summary, open(os.path.join(HERE, "pre_audit_results.json"), "w"), indent=2)
print("\n" + "=" * 64)
print("VERDICT:")
print(f"  - Bias (E-penalty) artifacts: 0 across all tiers (E-penalty drop validated).")
maxcorr = max(s["predicted_correlations"] for s in summary["tiers"].values())
print(f"  - Max predicted omit-C correlations in any tier: {maxcorr} "
      f"(vs RT where omit-C drove the C90-184 artifact).")
if concerns:
    print("  - ACTION: conserved cysteines left free (consider adding to fixed set):")
    for t, c, lab in concerns:
        print(f"      {t}: C{c} ({lab})")
else:
    print("  - No conserved cysteine is left free in any tier -> omit-C is safe as-is.")
print("\nwrote pre_audit_results.json")
