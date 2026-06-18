#!/usr/bin/env python3
"""Phylogeny-aware (Henikoff-reweighted) conservation tiers for SpRY-Cas9.

The default constraints/ stage computes conservation as a RAW, unweighted column
frequency, while its DCA stage already downweights redundant clades (neff 114 /
3782). That inconsistency lets clade over-sampling inflate "conservation". This
script recomputes conservation with the SAME per-sequence Henikoff weights used
for DCA, maps to PDB author numbering, and emits reweighted tiers + a diff vs the
unweighted tiers. The shared constraints/ module (and its tests) are left intact;
build_tiers.py picks up weighted_conservation_tiers.json when present.
"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from constraints.msa_io import parse_msa, encode_msa, filter_columns, AA_GAP, N_ALPHA
from constraints.coevolution import sequence_weights
from constraints.mapping import parse_structure_chain, build_column_to_pdb

MSA = os.path.join(REPO, "constraints/data/cas9_msa.fasta.gz")
CIF = os.path.join(REPO, "inputs/8SRS.cif")
SPRY = [61, 1111, 1135, 1136, 1218, 1219, 1317, 1322, 1333, 1335, 1337]
THRESHOLDS = (0.90, 0.70, 0.50)

print("loading MSA + computing Henikoff weights (matches DCA stage) ...")
seqs = parse_msa(MSA)                                   # list of (id, aligned_seq)
msa = encode_msa(seqs)                                  # (N, L)
wt_aligned = seqs[0][1]                                 # row-0 aligned sequence (NOT the tuple)
msa_f, _ = filter_columns(msa, max_gap_frac=0.5)
w = sequence_weights(msa_f, threshold=0.8)             # per-sequence weights
neff = float(w.sum() / w.max())
print(f"  N={msa.shape[0]} seqs, L={msa.shape[1]} cols, neff={neff:.1f} (vs raw {msa.shape[0]})")

# weighted top_frac per column (gaps excluded from denominator), vectorised
N, L = msa.shape
weighted_top = np.zeros(L)
for c in range(L):
    col = msa[:, c]
    nongap = col != AA_GAP
    denom = w[nongap].sum()
    if denom <= 0:
        continue
    # weighted count per amino-acid state (exclude gap index)
    wc = np.array([w[(col == a)].sum() for a in range(N_ALPHA - 1)])
    weighted_top[c] = wc.max() / denom

# map columns -> PDB author resnums (gap-safe; same mapping the pipeline uses)
struct = parse_structure_chain(CIF, chain="A")
col_to_pdb, _ = build_column_to_pdb(wt_aligned, struct, spry_positions=SPRY, strict=False)

def tiers_from(top_arr):
    out = {}
    for t in THRESHOLDS:
        out[f"ge_{int(t*100)}"] = sorted({col_to_pdb[c] for c in col_to_pdb if top_arr[c] >= t})
    return out

weighted_tiers = tiers_from(weighted_top)

# unweighted reference (from the existing constraints JSON) for the diff
con = json.load(open(os.path.join(REPO, "outputs/constraints/8srs_constraints.json")))
unw = con["conservation_tiers"]

print("\n=== conservation tier sizes: unweighted -> weighted (phylogeny-aware) ===")
diff_summary = {}
for k in ("ge_90", "ge_70", "ge_50"):
    a, b = set(unw[k]), set(weighted_tiers[k])
    entered, left = sorted(b - a), sorted(a - b)
    diff_summary[k] = {"unweighted": len(a), "weighted": len(b),
                       "entered": entered, "left": left}
    print(f"  {k}: {len(a):>4} -> {len(b):<4}  (+{len(entered)} entered, -{len(left)} left)")

out = {"method": "henikoff_reweighted", "threshold": 0.8, "neff": neff,
       "conservation_tiers": weighted_tiers,
       "conservation_tier_sizes": {k: len(v) for k, v in weighted_tiers.items()},
       "diff_vs_unweighted": diff_summary}
json.dump(out, open(os.path.join(HERE, "weighted_conservation_tiers.json"), "w"), indent=2)

# sanity: catalytic + highly-functional residues should stay strongly conserved
CAT = {10, 762, 839, 840, 863, 983, 986}
still_cat = [p for p in CAT if p in set(weighted_tiers["ge_50"])]
print(f"\ncatalytic residues still >=50% conserved after reweighting: "
      f"{len(still_cat)}/7 {sorted(still_cat)}")
print("wrote weighted_conservation_tiers.json")
