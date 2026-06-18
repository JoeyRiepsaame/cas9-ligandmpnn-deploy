#!/usr/bin/env bash
# Idempotent LigandMPNN setup for this machine (macOS / numpy-2 / no wget).
# Clones LigandMPNN, fetches the checkpoints we use via curl (get_model_params.sh
# needs wget, absent on macOS), and applies two required patches discovered during
# the 8SRS smoke test:
#   1. openfold uses removed numpy aliases (np.int/np.bool/np.object) -> numpy 2.x
#   2. run.py's per-design backbone-PDB write calls prody setResnames, which shape-
#      mismatches on the 8SRS nucleic-acid complex (5364 vs 5361). We don't need the
#      backbone PDBs (seqs + log_probs are in the stats .pt), so that block is removed.
# Safe to re-run.
set -euo pipefail
DIR="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/LigandMPNN}"
BASE="https://files.ipd.uw.edu/pub/ligandmpnn"
CKPTS=(proteinmpnn_v_48_020 solublempnn_v_48_020
       ligandmpnn_v_32_005_25 ligandmpnn_v_32_010_25 ligandmpnn_v_32_020_25)

if [[ ! -d "$DIR" ]]; then
  echo "[setup] cloning LigandMPNN -> $DIR"
  git clone -q https://github.com/dauparas/LigandMPNN.git "$DIR"
fi
mkdir -p "$DIR/model_params"
for ck in "${CKPTS[@]}"; do
  if [[ ! -s "$DIR/model_params/${ck}.pt" ]]; then
    echo "[setup] curl ${ck}.pt"
    curl -sf "$BASE/${ck}.pt" -o "$DIR/model_params/${ck}.pt"
  fi
done

# patch 1: numpy-2 deprecated aliases in openfold
perl -i -pe 's/np\.int(?![0-9a-zA-Z_])/int/g'      "$DIR/openfold/np/residue_constants.py"
perl -i -pe 's/np\.bool(?![0-9a-zA-Z_])/bool/g'    "$DIR/openfold/np/relax/utils.py"
perl -i -pe 's/np\.object(?![0-9a-zA-Z_])/object/g' "$DIR/openfold/data/templates.py"

# patch 2: remove the per-design backbone-PDB write block (setResnames mismatch)
if ! grep -q "backbone-PDB write disabled" "$DIR/run.py"; then
  python3 - "$DIR/run.py" <<'PY'
import sys
p=sys.argv[1]; src=open(p).read().splitlines(keepends=True); out=[]; i=0; n=len(src)
S="# write new sequences into PDB with backbone coordinates"; E="# write full PDB files"; done=False
while i<n:
    if S in src[i] and not done:
        ind=src[i][:len(src[i])-len(src[i].lstrip())]
        out.append(ind+"# [patched] backbone-PDB write disabled (seqs+log_probs in stats .pt;\n")
        out.append(ind+"# prody setResnames mismatches on 8SRS nucleic-acid complex)\n")
        i+=1
        while i<n and E not in src[i]: i+=1
        done=True; continue
    out.append(src[i]); i+=1
open(p,"w").write("".join(out)); print("[setup] run.py backbone block patched:", done)
PY
fi
python3 -c "import ast,sys; ast.parse(open('$DIR/run.py').read()); print('[setup] run.py OK')"
echo "[setup] LigandMPNN ready at $DIR"
