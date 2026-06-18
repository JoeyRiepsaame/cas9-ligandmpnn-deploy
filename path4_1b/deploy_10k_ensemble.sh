#!/usr/bin/env bash
# =============================================================================
# SpRY-Cas9 Path-4.1b — 10K single-constraint ENSEMBLE for LD robustness (gap C3)
# Separate from the per-tier SELECTION sweep (deploy_path4_1b.sh). Mirrors the RT
# campaign's 10K run: one constraint set (the operating tier, default T_balanced),
# big mixed-model ensemble, so design-ensemble LD/MI (ld_analysis.py) is stable.
# RT lesson: 60-design LD over-called couplings; 10K separates real from noise.
#
# Model mix (≈10,000 designs): SolubleMPNN 4000 + ProteinMPNN 3000 + LigandMPNN 3000
# Constraints: --omit_AA C, NO --bias_AA (E-penalty dropped), catalytic 7 fixed-WT.
# After this: extract_seqs -> ld_analysis.py --designs <ens> --map wt_resnum_to_index.json
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${ROOT:-$(dirname "$HERE")}"
LIGANDMPNN_DIR="${LIGANDMPNN_DIR:-$ROOT/LigandMPNN}"
INPUT_CIF="${INPUT_CIF:-$ROOT/inputs/8SRS.cif}"
TIER="${TIER:-T_balanced}"
FIXED="$(cat "$HERE/tiers/${TIER}_fixed.txt")"
OUT_DIR="${OUT_DIR:-$HERE/outputs/ensemble10k_${TIER}}"
DESIGN_CHAIN="A"

bash "$HERE/setup_ligandmpnn.sh" "$LIGANDMPNN_DIR"

# preflight: tier present + catalytic locked (reuse the same checks)
python3 - "$INPUT_CIF" "$HERE/tiers/${TIER}_fixed.txt" <<'PYEOF'
import sys
cif, fx = sys.argv[1], sys.argv[2]
chainsA=set(); in_atom=False; cols=[]; idx={}
for line in open(cif):
    if line.startswith('_atom_site.'): cols.append(line.strip().split('.')[1]); in_atom=True; continue
    if in_atom and line.startswith('ATOM'):
        if not idx: idx={c:i for i,c in enumerate(cols)}
        p=line.split()
        if p[idx['auth_asym_id']]=='A':
            try: chainsA.add(int(p[idx['auth_seq_id']]))
            except: pass
    elif in_atom and line.startswith('#'): in_atom=False
nums={int(t[1:]) for t in open(fx).read().split()}
assert not [n for n in nums if n not in chainsA], "fixed pos absent from chain A"
assert {10,762,839,840,863,983,986}<=nums, "catalytic not locked"
print(f"[ok] ensemble preflight: {len(nums)} fixed, catalytic locked, chain A {len(chainsA)} res")
PYEOF

gen () {  # gen <model_type> <ckpt_flag> <ckpt> <atom> <batch> <nbatch> <tag>
  local mtype="$1" cflag="$2" ckpt="$3" atom="$4" bs="$5" nb="$6" tag="$7"
  for temp in 0.1 0.3; do
    for seed in 1 2 3 4 5; do
      local out="$OUT_DIR/${tag}_T${temp}_s${seed}"; mkdir -p "$out"
      echo "[$TIER/10k] $mtype $tag T=$temp s=$seed  ${bs}x${nb}"
      python3 "$LIGANDMPNN_DIR/run.py" --model_type "$mtype" \
        $cflag "$LIGANDMPNN_DIR/model_params/${ckpt}.pt" \
        --pdb_path "$INPUT_CIF" --out_folder "$out" \
        --chains_to_design "$DESIGN_CHAIN" --fixed_residues "$FIXED" \
        --ligand_mpnn_use_atom_context "$atom" --omit_AA "C" \
        --temperature "$temp" --seed "$seed" \
        --batch_size "$bs" --number_of_batches "$nb" --save_stats 1
    done
  done
}

mkdir -p "$OUT_DIR"
# 2 temps x 5 seeds = 10 runs per line below; batch*nbatch sized to hit the targets:
gen soluble_mpnn --checkpoint_soluble_mpnn solublempnn_v_48_020 0 40 1 sol   # 10*40   = 400  ->x? see note
gen protein_mpnn --checkpoint_protein_mpnn proteinmpnn_v_48_020 0 30 1 prot  # 10*30   = 300
for ck in ligandmpnn_v_32_005_25 ligandmpnn_v_32_010_25 ligandmpnn_v_32_020_25; do
  gen ligand_mpnn --checkpoint_ligand_mpnn "$ck" 1 10 1 "lig_${ck##*_32_}"   # 3*10*10 = 300
done
# NOTE: per-line totals above are per the 10x(batch) loop. To reach ~10k scale up
# batch sizes (sol 400, prot 300, lig 100) for the production run; defaults here are
# a fast smoke profile (~1000) — set BATCH envs or edit before the real run.

echo ""
echo "DONE ensemble -> $OUT_DIR"
echo "Next: extract_seqs -> python3 ld_analysis.py --designs <ensemble.fasta|meta.json> \\"
echo "      --map wt_resnum_to_index.json --out-ld ld_matrix.json --out-score ld_scoring_table.json"
