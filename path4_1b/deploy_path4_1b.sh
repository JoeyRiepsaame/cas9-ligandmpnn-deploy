#!/usr/bin/env bash
# =============================================================================
# SpRY-Cas9 (8SRS) — Path-4.1b LigandMPNN generation (per-domain tiers)
# Ports the MMLV-RT Path-4.1b methodology onto SpRY-Cas9. Differences vs the
# original deploy_ligandmpnn_8srs.sh (which it supersedes):
#   1. NO  --bias_AA "E:-1.0"  -> the E-penalty is dropped; composition bias is
#      removed POST-HOC via Dutton correction (3 independent reasons, RT-4.1b).
#   2. MIX of model types (soluble + protein + ligand) -> RT-4.1b ESM ranking
#      showed Protein/Ligand designs top the naturalness axis, not only Soluble.
#   3. Per-DOMAIN tiers from path4_1b/tiers/ (PROTECT PI+HNH, REDESIGN REC).
#   4. All 7 catalytic residues are fixed-WT (active track). Nickase(H840A) and
#      dCas9(D10A+H840A) are derived POST-HOC as point mutations, not re-runs.
# Ranking here writes MPNN global_score only as ONE axis; final selection is
# windowed-ESM + Dutton-corrected MPNN + Pareto (see esm_score_windowed.py).
# =============================================================================
set -euo pipefail

# ---------- paths ----------
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${ROOT:-$(dirname "$HERE")}"                  # repo root (has inputs/, LigandMPNN/)
LIGANDMPNN_DIR="${LIGANDMPNN_DIR:-$ROOT/LigandMPNN}"
INPUT_CIF="${INPUT_CIF:-$ROOT/inputs/8SRS.cif}"
FIXED_DIR="${FIXED_DIR:-$HERE/tiers}"
OUT_DIR="${OUT_DIR:-$HERE/outputs}"                 # bulk designs (mirror to gdrive LigandMPNN)
DESIGN_CHAIN="A"
TIERS=(T_aggressive T_balanced T_safe)

# sampling per tier (~510/tier mirroring RT-4.1b selection sweep):
#   soluble : 2 temp x 3 seed x 20 = 120
#   protein : 2 temp x 3 seed x 20 = 120
#   ligand  : 3 ckpt x 2 temp x 3 seed x 15 = 270    -> 510 / tier, x3 tiers = 1530
TEMPS=(0.1 0.3)
SEEDS=(1 2 3)
LIGAND_CKPTS=("ligandmpnn_v_32_005_25" "ligandmpnn_v_32_010_25" "ligandmpnn_v_32_020_25")
PROTEIN_CKPT="proteinmpnn_v_48_020"
SOLUBLE_CKPT="solublempnn_v_48_020"
N_SOLUBLE=20; N_PROTEIN=20; N_LIGAND=15

# ---------- 0. setup (clone + curl checkpoints + numpy-2/run.py patches; idempotent) ----------
bash "$HERE/setup_ligandmpnn.sh" "$LIGANDMPNN_DIR"

# ---------- 1. PRE-FLIGHT (fail loudly) ----------
echo "=================== PRE-FLIGHT ==================="
python3 - "$INPUT_CIF" "$FIXED_DIR" "${TIERS[@]}" <<'PYEOF'
import sys, os
cif, fixed_dir = sys.argv[1], sys.argv[2]
tiers = sys.argv[3:]
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
assert {c for c in chains if c!='A'}, "FATAL: no nucleic-acid context -> atom_context no-op"
assert 'MG' in het, "WARN: catalytic Mg2+ not found"
print(f"[ok] chains={sorted(chains)}  HETATM={sorted(het)}")
print(f"[ok] chain A modeled residues: {len(chainsA)} (range {min(chainsA)}-{max(chainsA)})")
# every fixed position present; and the 7 catalytic must be fixed in every tier
CATALYTIC={10,762,839,840,863,983,986}
for tier in tiers:
    toks=open(os.path.join(fixed_dir,f"{tier}_fixed.txt")).read().split()
    nums={int(t[1:]) for t in toks}
    missing=[n for n in nums if n not in chainsA]
    assert not missing, f"FATAL: {tier} positions absent from chain A: {missing[:10]}"
    assert CATALYTIC<=nums, f"FATAL: {tier} missing catalytic {sorted(CATALYTIC-nums)}"
    print(f"[ok] {tier}: {len(nums)} fixed, all in chain A, 7/7 catalytic locked")
print("=== PRE-FLIGHT PASSED ===")
PYEOF

# ---------- 2. generation helpers ----------
# common flags: design chain A only, atom context, OMIT C, NO E-bias, save score+probs
gen () {  # gen <tier> <model_type> <ckpt_flag> <ckpt> <atom_ctx 0|1> <nseq> <tag>
  local tier="$1" mtype="$2" cflag="$3" ckpt="$4" atom="$5" nseq="$6" tag="$7"
  local fixed; fixed="$(cat "$FIXED_DIR/${tier}_fixed.txt")"
  for temp in "${TEMPS[@]}"; do
    for seed in "${SEEDS[@]}"; do
      local out="$OUT_DIR/$tier/${tag}_T${temp}_s${seed}"; mkdir -p "$out"
      echo "[$tier] $mtype $tag T=$temp seed=$seed n=$nseq -> $out"
      python3 "$LIGANDMPNN_DIR/run.py" \
        --model_type "$mtype" \
        $cflag "$LIGANDMPNN_DIR/model_params/${ckpt}.pt" \
        --pdb_path "$INPUT_CIF" --out_folder "$out" \
        --chains_to_design "$DESIGN_CHAIN" --fixed_residues "$fixed" \
        --ligand_mpnn_use_atom_context "$atom" \
        --omit_AA "C" \
        --temperature "$temp" --seed "$seed" \
        --batch_size "$nseq" --number_of_batches 1 \
        --save_stats 1
    done
  done
}

run_tier () {
  local tier="$1"
  echo ""; echo "############### $tier ###############"
  gen "$tier" soluble_mpnn --checkpoint_soluble_mpnn "$SOLUBLE_CKPT" 0 "$N_SOLUBLE" "sol"
  gen "$tier" protein_mpnn --checkpoint_protein_mpnn "$PROTEIN_CKPT" 0 "$N_PROTEIN" "prot"
  for ck in "${LIGAND_CKPTS[@]}"; do
    gen "$tier" ligand_mpnn --checkpoint_ligand_mpnn "$ck" 1 "$N_LIGAND" "lig_${ck##*_32_}"
  done
}

mkdir -p "$OUT_DIR"
for t in "${TIERS[@]}"; do run_tier "$t"; done

# ---------- 3. flat MPNN-score table (one axis; ESM+Dutton+Pareto come later) ----------
echo ""; echo "=== collecting MPNN global_score (one ranking axis) ==="
python3 - "$OUT_DIR" "${TIERS[@]}" <<'PYEOF'
import sys, glob, os, re, csv
root=sys.argv[1]; tiers=set(sys.argv[2:])
rows=[]
for fa in glob.glob(os.path.join(root,'**','seqs','*.fa'), recursive=True):
    parts=os.path.relpath(fa,root).split(os.sep); tier=parts[0]
    for line in open(fa):
        if line.startswith('>'):
            m=re.search(r'overall_confidence=([0-9.]+)', line)
            if m: rows.append((tier, float(m.group(1)), os.path.relpath(fa,root), line.strip()))
rows.sort(key=lambda r:(r[0], r[1]))
with open(os.path.join(root,'ranked_designs.csv'),'w',newline='') as fh:
    w=csv.writer(fh); w.writerow(['tier','global_score','source','header']); w.writerows(rows)
from collections import Counter
print("designs per tier:", dict(Counter(r[0] for r in rows)), " total:", len(rows))
print("wrote", os.path.join(root,'ranked_designs.csv'))
PYEOF

echo ""
echo "DONE. Next: extract_seqs -> esm_score_windowed.py -> dutton_correct -> pareto -> AF3."
echo "Nickase/dCas9 variants are post-hoc point mutations on the active-track winners."
