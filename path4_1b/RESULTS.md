# SpRY-Cas9 Path-4.1b — Generation + Scoring Results

**Date:** 2026-06-18 · local MPS run · target: SpRY-Cas9 (8SRS chain A, 1341 aa)

## Pipeline executed
generation (1,530 designs) → extract/dedup → Dutton-corrected MPNN → windowed ESM-2
→ ESM⟂Dutton Pareto + %identity gate → 17-design shortlist → design-ensemble LD.

## Generation
- **1,530 designs** = 510 / tier × 3 tiers (Soluble+Protein+Ligand mix), `--omit_AA C`, no E-penalty, catalytic fixed-WT.
- **100% unique**; omit-C perfect (cysteine only at fixed C574 in every design).

## Per-tier frontier (the ESM ⟂ MPNN trade-off, RT-consistent)
| Tier | mean ESM (naturalness) | mean Dutton-MPNN | mean %id to WT | patent <70% |
|---|---|---|---|---|
| T_aggressive | −0.371 (lowest) | **2.735 (highest)** | **57.0** | ✓ clears easily |
| T_balanced | −0.333 | 2.610 | 68.2 | ✓ (borderline) |
| T_safe | **−0.318 (highest)** | 2.520 | 75.4 | ✗ fails gate |

WT ESM = −0.307 (all designs below WT naturalness, expected at this divergence).
The axes are **anti-correlated**: most-divergent tier wins MPNN, most-fixed tier wins
naturalness — exactly why selection uses the Pareto frontier, not either axis alone.

## Cross-model fairness fix (per-model-z)
The campaign generated all 3 model types (Ligand 810 / Protein 360 / Soluble 360) — a real
competition (unlike the RT sweep, which was LigandMPNN-only). BUT the Dutton/MPNN score is
**not comparable across model types**: LigandMPNN conditions on the DNA/RNA atom context, so it
assigns systematically higher self-confidence (mean mpnn_corr Ligand 2.656 > Protein 2.600 >
Soluble 2.565) — an information advantage, not a quality one. On the model-agnostic **ESM** axis
the three are **near-tied** (−0.339 / −0.341 / −0.343). The naive Pareto therefore swept 17/17
LigandMPNN. **Fix:** z-normalise the Dutton score WITHIN each model type before the cross-model
Pareto, so Soluble/Protein compete fairly. (`pareto_rank.py`, per-model-z.)

## Synthesis shortlist — 20 designs (per-model-z Pareto ∩ <70% identity)
- **By model: 13 Ligand + 5 Protein + 2 Soluble** (Ligand's 65% ≈ its 53% pool share — fair,
  vs the 100% sweep before the fix). The fix surfaced 7 Protein/Soluble designs the biased axis hid.
- **By tier: 11 T_aggressive (~57–61% id) + 9 T_balanced (~69%)** — bold-vs-safe flavors retained.
- File: `synthesis_shortlist.json`. **Re-audited:** all 20 byte-match source, len 1341,
  7/7 catalytic WT, 11/11 SpRY preserved, omit-C clean.
- AF3 batch updated to 25 jobs (20 designs + 5 calibrators).

## Design-ensemble LD (artifact scan)
- 729 variable positions; mean MI-APC ≈ 0.0002 (ensemble largely independent → healthy).
- **No omit-C / bias artifact** — no top coupling involves C80 (contrast RT's C90–184); E-penalty dropped.
- Mild coupling hub at A122 (pairs with 519/721/1349) — flag, not a blocker.

## Interpretation / next steps
- Two operating flavors: **T_aggressive (~60% id)** for maximal novelty (function-riskier) vs
  **T_balanced (~69%)** for safer function (borderline novelty). The <70% multidomain-RT frontier
  was unproven in the lit review, so both warrant a wet-lab pilot.
- **AF3 triage (Phase 6):** fold/binding sanity gate only (with the calibrator ladder) — AF3 cannot
  rank function or test PAM relaxation.
- **Wet-lab (Phase 8):** the real test — cleavage on NGG + NAA/NAC/NAG PAMs.
- Nickase (H840A) / dCas9 (D10A+H840A) variants are post-hoc point mutations on the chosen winners.
