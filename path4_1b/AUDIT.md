# SpRY-Cas9 Path-4.1b — Extensive Validation Audit

**Date:** 2026-06-19 · **95 independent checks across 3 validators, 0 failures.**
All re-derive ground truth (structure, MSA, .pt files, AF3 outputs) rather than trusting derived JSONs.

| Validator | Checks | Scope | Result |
|---|---|---|---|
| `validate_results.py` | 46 | foundation: numbering, tiers, conservation, constants | ✅ 0 fail |
| `validate_pipeline.py` | 38 | generation, scoring, per-model-z, shortlists, pilot, AF3 batch | ✅ 0 fail |
| `validate_af3.py` | 11 | AF3 round-trip + metric re-extraction + conclusion | ✅ 0 fail |

## Foundation (46) — re-derived vs canonical SpCas9 / MSA / structure
- Numbering == canonical SpCas9 (Q99ZW2), 0 offsets across 1341 residues; SpRY/catalytic exact.
- Constants (11 SpRY + 7 catalytic) identical across all artifacts; WT byte-matches structure.
- Unweighted conservation reproduces validated 490/741/1048; phylogeny-aware tiers nested + catalytic-safe.
- Tiers: always-fixed(95) ⊆ each; nesting Tagg⊂Tbal⊂Tsafe; deterministic regeneration.

## Pipeline + data (38) — independent decode of all 1530 designs
- **All 1530 designs decoded straight from the 90 stats .pt**: byte-match `unique_meta`;
  **every design** len-1341 / omit-C / 7-catalytic / 11-SpRY (0 violations).
- Dutton background + raw MPNN reproduced from .pt; ESM covers 1530, no NaN.
- **Per-model-z fix verified** (z_mpnn == per-model normalised, 0 mismatches); Ligand atom-context
  bias confirmed (2.656 > 2.600 > 2.565).
- Shortlists: per-model-z all on Pareto frontier, all <70 %, span all 3 models; overlap=4 reproduced.
- Pilot: consensus ⊆ pilot, byte-match; AF3 batch layout/count/byte-match.

## AF3 results (11) — round-trip + headline conclusion
- **Every one of the 22 folded proteins BYTE-MATCHES the source sequence** — AF3 folded exactly
  what we sent (designs vs pilot_set, calibrators vs af3_calibrators).
- 6-chain 8SRS layout in every job; Cas9:sgRNA + R-loop re-extracted from summary JSONs match the table.
- **scram_rec non-discrimination reproduced independently:** local pLDDT @200–260 = scram_rec 43.7,
  native dCas9 88.6, real designs 41.8/40.7/47.5 → AF3 cannot tell a heavy REC redesign from
  random nonsense (the calibrator ladder's key finding, validated).

## Bugs caught + fixed earlier this campaign (all now regression-covered)
tuple-vs-seq mapping · int8 overflow at scale · `open()[]`-index · openfold numpy-2 · `--save_stats`
flag · setResnames crash · cross-model MPNN-score bias (per-model-z) · AF3 SpRY-R61 calibrator slip.

## Verdict
The pipeline and data are **internally consistent and match ground truth end-to-end**. The 18-design
pilot is sequence-audited and fairly ranked; AF3 is confirmed non-discriminating (a documented,
validated limitation, not a pipeline error). Wet-lab PAM panel is the sole functional arbiter.
