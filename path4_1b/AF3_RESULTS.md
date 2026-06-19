# SpRY-Cas9 Path-4.1b — AF3 Triage Results (complete)

**Date:** 2026-06-19 · folded **20/23 jobs: 17 designs + 3 calibrators** (bind_dead, dcas9, ncas9_h840a).
Still pending: **wt_spry** (positive) and **scram_rec** (fold-negative). Table: `af3_interface_table.csv`.

## Calibrators anchor the WT-fold band
All 3 calibrators fold/bind near-identically: **pLDDT ~89.7, Cas9:sgRNA ~0.86, Cas9:DNA ~0.95,
cat_min ~91.1**. dCas9 (D10A+H840A) folds/binds exactly like WT → it serves as the **WT-fold
positive anchor**. This is the "good" reference band.

## R-loop iPTM = artifact (now confirmed by dCas9)
sgRNA:DNA iPTM = 0.03 for **dCas9 too** — dCas9 binds perfectly, yet scores the same low R-loop as
every design. Cause: 98 nt sgRNA dilutes the whole-chain iPTM + AF3 keeps the DNA duplexed
(DNA:DNA 0.24–0.33 > sgRNA:DNA 0.01–0.03). **Do not rank on R-loop.** R-loop/PAM → wet lab.

## All 17 designs PASS the fold + binding gate (but below the WT ceiling)
Every design: iPTM 0.90–0.92, Cas9:DNA 0.93–0.94, catalytic site ordered. They sit ~5 pLDDT
below the calibrator ceiling (84–86 vs 89.7) and lower cat_min — expected for 30–43 % redesign,
none catastrophic.

## Ranking vs the WT-fold band (pLDDT + Cas9:gRNA + catalytic ordering)
| Tier | designs | note |
|---|---|---|
| **TOP** | **d00806\*, d00900\*, d00941** | closest to WT-fold on all 3 metrics; **both consensus designs are TOP** |
| mid | d00107, d00824, d00918, d01011, d00526, d00523, d00537, d00097, d00018 | solid fold; d00537 highest pLDDT (86.1); d00107 highest gRNA (0.88)+cat |
| **weak — deprioritise** | d00272 (pLDDT 80.7), d00542\* (gRNA 0.75), d00134 (0.75), d00102 (gRNA 0.75, cat 75.4), d00101 (gRNA 0.70) | most are the most-divergent bold/T_aggressive picks |

\*consensus. Reassuring: the two designs both scoring axes agreed on (d00806, d00900) fold best.

## Remaining gap + recommendation
- **Only `scram_rec` (fold-negative) is still critical** — it sets the floor; without it we can't be
  sure AF3 discriminates fold quality (vs scoring everything ~84). dCas9 already gives the positive.
  → fold `scram_rec` (1 job) to close the gate. (wt_spry optional — dCas9 ≈ WT-fold.)
- **Carry forward:** the TOP set (d00806, d00900, d00941) + a couple of mids (d00107 strong gRNA,
  d00537 high pLDDT) + 1–2 bold (d00097) for the divergence test. Deprioritise d00101/d00102/d00272.
- **Wet-lab PAM panel = the real arbiter** (AF3 can't test catalysis or PAM relaxation).
