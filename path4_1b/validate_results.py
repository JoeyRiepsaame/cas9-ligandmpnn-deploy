#!/usr/bin/env python3
"""THOROUGH end-to-end validation of all path4_1b results.

Independent re-verification: instead of trusting each script's own output, this
re-derives the key facts from authoritative sources (8SRS structure, the MSA,
constraints/pipeline.py) and cross-checks every artifact against ground truth AND
against each other. Then it invokes the component self-tests. Exit non-zero on any
failure. Two bugs were already caught this session (a tuple-vs-seq mapping bug and
an int8 overflow), so this is deliberately paranoid.

Run:  python validate_results.py
"""
import json, os, re, sys, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO); sys.path.insert(0, HERE)

FAIL, NCHECK = [], 0
def ck(cond, msg):
    global NCHECK; NCHECK += 1
    print(("  PASS  " if cond else "  FAIL  ") + msg)
    if not cond: FAIL.append(msg)
def sec(t): print("\n" + "=" * 72 + f"\n{t}\n" + "=" * 72)

def load(p): return json.load(open(os.path.join(HERE, p)))
def fixed_set(name):
    return {int(t[1:]) for t in open(os.path.join(HERE, "tiers", name)).read().split()}
def fasta(p):
    return "".join(l.strip() for l in open(p) if not l.startswith(">"))

# ---- authoritative ground truth (re-derived, not trusted from our own JSON) ----
from constraints.mapping import parse_structure_chain, build_column_to_pdb
from constraints.conservation import compute_conservation, conserved_positions
# SpRY from the canonical pipeline module; catalytic from the original repo tier1.
from constraints.pipeline import SPRY_POSITIONS
SPRY = set(SPRY_POSITIONS)
CATALYTIC = {10, 762, 839, 840, 863, 983, 986}
struct = parse_structure_chain(os.path.join(REPO, "inputs/8SRS.cif"), chain="A")
resnums = sorted(struct)
wt_struct = "".join(struct[r] for r in resnums)

# =====================================================================
sec("1. GROUND-TRUTH CONSTANTS — consistency across every artifact")
dm = load("domain_map.json")
ck(set(s["pos"] for s in dm["spry_positions"]) == SPRY, "domain_map SpRY == pipeline SPRY_POSITIONS")
ck(set(c["pos"] for c in dm["catalytic"]) == CATALYTIC, "domain_map catalytic == {10,762,839,840,863,983,986}")
con = load(os.path.join(REPO, "outputs/constraints/8srs_constraints.json"))
ck(set(con["spry_positions"]) == SPRY, "constraints JSON spry_positions == SPRY")
pa = load("pre_audit_results.json")
# every tier file must contain SpRY and catalytic
for t in ("T_aggressive", "T_balanced", "T_safe"):
    fx = fixed_set(f"{t}_fixed.txt")
    ck(SPRY <= fx and CATALYTIC <= fx, f"{t}: contains all 11 SpRY + 7 catalytic")
# deploy scripts hardcode the catalytic preflight set — check it matches
for sh in ("deploy_path4_1b.sh", "deploy_10k_ensemble.sh"):
    txt = open(os.path.join(HERE, sh)).read()
    m = re.search(r"\{10,762,839,840,863,983,986\}", txt)
    ck(m is not None, f"{sh}: catalytic preflight set matches ground truth")

# =====================================================================
sec("2. STRUCTURE & WT REFERENCE")
ck(len(wt_struct) == 1341, f"8SRS chain A modeled residues == 1341 (got {len(wt_struct)})")
ck((resnums[0], resnums[-1]) == (3, 1366), f"author range == 3-1366 (got {resnums[0]}-{resnums[-1]})")
wt_file = fasta(os.path.join(HERE, "spry_cas9_wt.fasta"))
ck(wt_file == wt_struct, "spry_cas9_wt.fasta byte-matches structure-derived sequence")
r2i = {int(k): v for k, v in load("wt_resnum_to_index.json").items()}
ck(len(r2i) == 1341 and set(r2i) == set(resnums), "wt_resnum_to_index covers all 1341 modeled resnums")
ck(all(r2i[r] == i for i, r in enumerate(resnums)), "index map is monotonic with sorted author resnums")
ck(all(struct[p] == lab[0] for p, lab in [(c["pos"], c["label"]) for c in dm["catalytic"]]),
   "catalytic residues carry expected AA in structure")
known_spry = {61:"R",1111:"R",1135:"L",1136:"W",1218:"K",1219:"Q",1317:"R",1322:"R",1333:"P",1335:"Q",1337:"R"}
ck(all(struct[p] == aa for p, aa in known_spry.items()), "SpRY positions carry expected substitution residue")

# =====================================================================
sec("3. CONSERVATION — unweighted reproducibility + phylogeny-aware re-rank")
# unweighted must reproduce the validated 490/741/1048 from the constraints stage
seqs_msa = None
cons = compute_conservation(os.path.join(REPO, "constraints/data/cas9_msa.fasta.gz"))
from constraints.msa_io import parse_msa
wt_aligned = parse_msa(os.path.join(REPO, "constraints/data/cas9_msa.fasta.gz"))[0][1]
c2p, _ = build_column_to_pdb(wt_aligned, struct, spry_positions=list(SPRY), strict=False)
for thr, exp in [(0.90, 490), (0.70, 741), (0.50, 1048)]:
    got = len(conserved_positions(cons, c2p, thr))
    ck(got == exp, f"unweighted conservation >= {int(thr*100)}% == {exp} (got {got})")
# weighted tiers: catalytic retained, sizes self-consistent, re-rank diff correct
wt_tiers = load("weighted_conservation_tiers.json")
wts = {k: set(v) for k, v in wt_tiers["conservation_tiers"].items()}
ck(CATALYTIC <= wts["ge_50"], "phylogeny-aware: all 7 catalytic still >=50% conserved")
ck(wts["ge_90"] <= wts["ge_70"] <= wts["ge_50"], "weighted tiers are nested (ge90 subset ge70 subset ge50)")
ck(all(wt_tiers["conservation_tier_sizes"][k] == len(wts[k]) for k in wts),
   "weighted tier sizes match their lists")
# verify diff_vs_unweighted internally consistent
unw = {k: set(con["conservation_tiers"][k]) for k in wts}
for k in wts:
    d = wt_tiers["diff_vs_unweighted"][k]
    ck(set(d["entered"]) == wts[k] - unw[k] and set(d["left"]) == unw[k] - wts[k],
       f"{k}: recorded entered/left sets match recomputed weighted-vs-unweighted diff")

# =====================================================================
sec("4. TIERS — validity, always-fixed, nesting, determinism")
modeled = set(resnums)
contacts = fixed_set("../fixed_residues/tier1_fixed.txt".replace("tiers/", "")) if False else \
           ({int(t[1:]) for t in open(os.path.join(REPO,"fixed_residues/tier1_fixed.txt")).read().split()} - SPRY - CATALYTIC)
ALWAYS = SPRY | CATALYTIC | contacts | set(con["dca_top"])
ck(len(contacts) == 35, f"contacts derived = tier1 - SpRY - catalytic == 35 (got {len(contacts)})")
TA, TB, TS = (fixed_set(f"{t}_fixed.txt") for t in ("T_aggressive","T_balanced","T_safe"))
for name, fx in [("T_aggressive",TA),("T_balanced",TB),("T_safe",TS)]:
    ck(fx <= modeled, f"{name}: all fixed positions exist in chain A")
    ck(ALWAYS <= fx, f"{name}: always-fixed (SpRY+cat+contacts+DCA={len(ALWAYS)}) subset of tier")
ck(TA <= TB <= TS, "tier nesting T_aggressive subset T_balanced subset T_safe")
summ = load("tiers_summary.json")
ck(all(summ["tiers"][n]["n_fixed"] == len(fx)
       for n, fx in [("T_aggressive",TA),("T_balanced",TB),("T_safe",TS)]),
   "tiers_summary n_fixed matches the actual tier files")
# determinism: regenerate tiers and compare set-equality (does not trust on-disk copy)
snap = {n: set(fx) for n, fx in [("T_aggressive",TA),("T_balanced",TB),("T_safe",TS)]}
r = subprocess.run([sys.executable, "build_tiers.py"], cwd=HERE, capture_output=True, text=True)
ck(r.returncode == 0, "build_tiers.py re-runs cleanly")
ck("phylogeny-aware" in r.stdout, "build_tiers consumes the phylogeny-aware (weighted) conservation")
regen = {n: fixed_set(f"{n}_fixed.txt") for n in ("T_aggressive","T_balanced","T_safe")}
ck(regen == snap, "tiers are deterministic (regeneration identical to committed files)")

# =====================================================================
sec("5. COMPONENT SELF-TESTS / AUDITS (subprocess)")
def run(label, args, expect="PASS"):
    r = subprocess.run([sys.executable] + args, cwd=HERE, capture_output=True, text=True)
    ok = r.returncode == 0 and (expect in (r.stdout + r.stderr) if expect else True)
    ck(ok, f"{label} (exit {r.returncode}{', '+expect+' found' if ok and expect else ''})")
    return r
run("audit_numbering.py", ["audit_numbering.py"], "AUDIT PASSED")
run("esm_score_windowed --selftest", ["esm_score_windowed.py", "--selftest"], "SELFTEST PASS")
run("ld_analysis --selftest", ["ld_analysis.py", "--selftest"], "SELFTEST PASS")
rp = subprocess.run([sys.executable, "pre_audit.py"], cwd=HERE, capture_output=True, text=True)
ck(rp.returncode == 0 and "omit-C is safe as-is" in rp.stdout, "pre_audit re-runs and verdict = omit-C safe")

# =====================================================================
sec("6. DEPLOY SCRIPTS — syntax + no E-penalty (Path-4.1b decision)")
for sh in ("deploy_path4_1b.sh", "deploy_10k_ensemble.sh"):
    r = subprocess.run(["bash", "-n", os.path.join(HERE, sh)], capture_output=True, text=True)
    ck(r.returncode == 0, f"{sh}: bash syntax OK")
    txt = open(os.path.join(HERE, sh)).read()
    # strip full-line comments so we test ACTUAL flag usage, not documentation of it
    code = "\n".join(ln for ln in txt.splitlines() if not ln.lstrip().startswith("#"))
    ck("--omit_AA" in code and '"C"' in code, f"{sh}: omit_AA C present")
    ck("--bias_AA" not in code, f"{sh}: NO --bias_AA flag in any executable line (E-penalty dropped)")

# =====================================================================
print("\n" + "=" * 72)
print(f"VALIDATION: {NCHECK} checks, {len(FAIL)} failures")
if FAIL:
    for m in FAIL: print("   FAIL: " + m)
    sys.exit(1)
print("ALL RESULTS VALIDATED — artifacts internally consistent and match ground truth.")
