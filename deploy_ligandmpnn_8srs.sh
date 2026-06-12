#!/usr/bin/env bash
# =============================================================================
# SpRY-Cas9 (8SRS) — 3-tier LigandMPNN deploy for RunPod
# Target: chain A (Cas9, 1341 modeled residues, auth numbering 3-1366)
# Context: chain B (gRNA), C/c/D (target/non-target DNA), MG2+  -> atom context
# Evidence basis: 12-paper audit (2026-06-12). Conservation backbone + 35 DCA
#                 positions + active-site contacts + PTM + SpRY. See audit .md.
#
# Tiers (audited, nested T1 subset T2 subset T3):
#   T1 =   53 fixed ( 4.0%)  active-site contacts + SpRY only  -> expression control
#   T2 =  692 fixed (51.6%)  consv>=90% + 35 DCA + contacts + PTM + SpRY -> most novel
#   T3 = 1120 fixed (83.5%)  consv>=50% + DCA + contacts + PTM + SpRY    -> safest (=Tao PE8 range)
#
# Sampling per tier: 3 checkpoints x 2 temps x 2 seeds x 10 seqs = 120 designs
# Ranking: MPNN log-likelihood (--save_score 1), best single predictor (Johnson 2025)
# =============================================================================
set -euo pipefail

# ---------- paths ----------
ROOT="${ROOT:-$PWD}"
LIGANDMPNN_DIR="${LIGANDMPNN_DIR:-$ROOT/LigandMPNN}"
INPUT_CIF="${INPUT_CIF:-$ROOT/inputs/8SRS.cif}"
FIXED_DIR="${FIXED_DIR:-$ROOT/fixed_residues}"
OUT_DIR="${OUT_DIR:-$ROOT/outputs}"
DESIGN_CHAIN="A"

# 3 ligand-MPNN checkpoints = 3 training noise levels (0.05 / 0.10 / 0.20 Å) for diversity
CKPTS=(
  "ligandmpnn_v_32_005_25"
  "ligandmpnn_v_32_010_25"
  "ligandmpnn_v_32_020_25"
)
TEMPS=(0.1 0.3)
SEEDS=(1 2)
SEQS_PER_RUN=10            # batch_size; 3 ckpt x 2 temp x 2 seed x 10 = 120 / tier

# ---------- 0. setup ----------
if [[ ! -d "$LIGANDMPNN_DIR" ]]; then
  echo "[setup] cloning LigandMPNN ..."
  git clone https://github.com/dauparas/LigandMPNN.git "$LIGANDMPNN_DIR"
  ( cd "$LIGANDMPNN_DIR" && bash get_model_params.sh "./model_params" )
  pip install -q torch numpy ProDy biopython ml-collections
fi

# ---------- 1. PRE-FLIGHT VALIDATION (fail loudly, not silently) ----------
echo "=================== PRE-FLIGHT ==================="
python3 - "$INPUT_CIF" "$FIXED_DIR" <<'PYEOF'
import sys
cif, fixed_dir = sys.argv[1], sys.argv[2]

# (a) structure exists & has the nucleic-acid/ligand context
chainsA=set(); chains=set(); het=set()
with open(cif) as f:
    in_atom=False; cols=[]; idx={}
    for line in f:
        if line.startswith('_atom_site.'):
            cols.append(line.strip().split('.')[1]); in_atom=True; continue
        if in_atom and (line.startswith('ATOM') or line.startswith('HETATM')):
            if not idx: idx={c:i for i,c in enumerate(cols)}
            p=line.split(); a=p[idx['auth_asym_id']]; chains.add(a)
            if line.startswith('ATOM') and a=='A':
                try: chainsA.add(int(p[idx['auth_seq_id']]))
                except: pass
            if line.startswith('HETATM'): het.add(p[idx['auth_comp_id']])
        elif in_atom and line.startswith('#'): in_atom=False

assert 'A' in chains, "FATAL: design chain A absent"
nuc = {c for c in chains if c!='A'}
assert nuc, "FATAL: no nucleic-acid context chains -> --ligand_mpnn_use_atom_context is a no-op"
assert 'MG' in het, "WARN: catalytic Mg2+ not found"
print(f"[ok] chains={sorted(chains)}  HETATM={sorted(het)}")
print(f"[ok] chain A modeled residues: {len(chainsA)} (range {min(chainsA)}-{max(chainsA)})")

# (b) every fixed position must exist in chain A author numbering
import os
for tier in ('tier1','tier2','tier3'):
    toks=open(os.path.join(fixed_dir,f"{tier}_fixed.txt")).read().split()
    nums=[int(t[1:]) for t in toks]
    missing=[n for n in nums if n not in chainsA]
    assert not missing, f"FATAL: {tier} positions absent from chain A: {missing[:10]}"
    print(f"[ok] {tier}: {len(nums)} fixed positions all present in chain A")
print("=== PRE-FLIGHT PASSED ===")
PYEOF

# ---------- 2. RUN DESIGNS ----------
run_tier () {
  local tier="$1"
  local fixed; fixed="$(cat "$FIXED_DIR/${tier}_fixed.txt")"
  echo ""
  echo "############### $tier ###############"
  for ckpt in "${CKPTS[@]}"; do
    for temp in "${TEMPS[@]}"; do
      for seed in "${SEEDS[@]}"; do
        local out="$OUT_DIR/$tier/${ckpt}_T${temp}_s${seed}"
        mkdir -p "$out"
        echo "[$tier] ckpt=$ckpt T=$temp seed=$seed -> $out"
        python3 "$LIGANDMPNN_DIR/run.py" \
          --model_type "ligand_mpnn" \
          --checkpoint_ligand_mpnn "$LIGANDMPNN_DIR/model_params/${ckpt}.pt" \
          --pdb_path "$INPUT_CIF" \
          --out_folder "$out" \
          --chains_to_design "$DESIGN_CHAIN" \
          --fixed_residues "$fixed" \
          --ligand_mpnn_use_atom_context 1 \
          --omit_AA "C" \
          --bias_AA "E:-1.0" \
          --temperature "$temp" \
          --seed "$seed" \
          --batch_size "$SEQS_PER_RUN" \
          --number_of_batches 1 \
          --save_score 1 \
          --save_probs 1
      done
    done
  done
}

mkdir -p "$OUT_DIR"
run_tier tier1
run_tier tier2
run_tier tier3

# ---------- 3. RANK BY MPNN SCORE ----------
echo ""
echo "=== ranking all designs by MPNN log-likelihood (global_score) ==="
python3 - "$OUT_DIR" <<'PYEOF'
import sys, glob, os, re
root=sys.argv[1]
rows=[]
for fa in glob.glob(os.path.join(root,'**','seqs','*.fa'), recursive=True):
    tier=fa.split(os.sep)[len(root.split(os.sep))]
    name=None
    for line in open(fa):
        if line.startswith('>'):
            m=re.search(r'global_score=([0-9.]+)', line)   # lower = better
            if m:
                rows.append((tier, float(m.group(1)), os.path.relpath(fa,root), line.strip()))
rows.sort(key=lambda r:(r[0], r[1]))
print(f"{'tier':6} {'score':>8}  source")
for t,s,src,_ in rows[:40]:
    print(f"{t:6} {s:8.4f}  {src}")
print(f"\nTotal designs scored: {len(rows)}")
# write a flat csv
import csv
with open(os.path.join(root,'ranked_designs.csv'),'w',newline='') as fh:
    w=csv.writer(fh); w.writerow(['tier','global_score','source','header'])
    w.writerows(rows)
print("Wrote", os.path.join(root,'ranked_designs.csv'))
PYEOF

echo ""
echo "DONE. Designs in $OUT_DIR/<tier>/...  Ranked: $OUT_DIR/ranked_designs.csv"
echo "Next: take top-scoring designs per tier -> AF3 structural triage -> order panel."
