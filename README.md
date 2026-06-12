# SpRY-Cas9 (8SRS) — 3-Tier LigandMPNN Deploy

**Date:** 2026-06-12 · Target: SpRY-Cas9 chain A (8SRS) · Run on: RunPod GPU

## What this is
Evidence-based 3-tier LigandMPNN redesign of SpRY-Cas9, using the 8SRS
`SpRY-Cas9:gRNA:DNA` complex so the gRNA + target/non-target DNA + Mg²⁺ act as
**atom context** (`--ligand_mpnn_use_atom_context 1`). Constraint strategy comes
from a 12-paper wet-lab audit (see `../SpRY-Cas9-MPNN-constraint-audit-...md`).

## Tiers (audited, nested T1 ⊂ T2 ⊂ T3)
| Tier | Fixed | % | Composition | Purpose |
|---|---|---|---|---|
| T1 | 53 | 4.0% | active-site contacts + SpRY | expression control (predicted no cut) |
| T2 | 692 | 51.6% | consv≥90% + 35 DCA + contacts + PTM + SpRY | most novel; tests DCA |
| T3 | 1120 | 83.5% | consv≥50% + DCA + contacts + PTM + SpRY | safest (= Tao PE8 64–84% range) |

## Files
```
inputs/8SRS.cif            # structure (mmCIF — preserves lowercase chain 'c')
fixed_residues/tierN_fixed.txt   # validated "A53 A61 ..." fixed-position strings
deploy_ligandmpnn_8srs.sh  # main script (clones LigandMPNN, preflight, run, rank)
```

## Run on RunPod
```bash
cd SpRY-Cas9-LigandMPNN-deploy
bash deploy_ligandmpnn_8srs.sh         # first run clones LigandMPNN + model params
```
Produces **120 designs/tier** (3 checkpoints × 2 temps[0.1,0.3] × 2 seeds × 10 seqs)
and `outputs/ranked_designs.csv` sorted by MPNN `global_score` (lower = better;
best single predictor of in-vitro activity per Johnson 2025).

## Built-in safety guards (preflight — script aborts if any fail)
- Confirms chain A present + nucleic-acid context chains exist (else atom-context is a no-op)
- Confirms catalytic Mg²⁺ present
- Asserts **every** fixed position exists in chain A author numbering (3–1366, 23 gaps) — no silent drops

## Verified pre-flight (2026-06-12)
- 8SRS chains: A(Cas9 1341 aa) · B(gRNA) · C/c/D(DNA) · MG²⁺ ✓
- All T1/T2/T3 positions present in chain A — 0 missing ✓
- bash syntax OK ✓

## Settings rationale
- `--omit_AA "C"` avoid spurious cysteines · `--bias_AA "E:-1.0"` mild Glu down-bias
- `--chains_to_design A` redesign Cas9 only; B/C/D/c/Mg = fixed context
- `--save_score 1 --save_probs 1` for MPNN-score ranking

## Next
top-scoring designs/tier → AF3 structural triage → synthesis panel.
⚠️ Per lab policy: any sequence leaving this pipeline must pass the file-based
sequence-verification audit before going to chat or gene synthesis.
