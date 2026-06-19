# SpRY-Cas9 Path-4.1b — AF3 Triage Results

**Date:** 2026-06-19 · folded 16/23 jobs: **14 designs + 2 calibrators** (bind_dead, ncas9_h840a).
Table: `af3_interface_table.csv`. Chain order [0]Cas9 [1]sgRNA [2-4]DNA [5]Mg.

## Headline — all 14 designs PASS the fold + substrate-binding sanity gate
- **Fold:** iPTM 0.90–0.92, pTM 0.82–0.87, pLDDT 80.7–85.7 → every design folds confidently.
- **Substrate binding:** Cas9:DNA 0.93–0.95 (all), Cas9:sgRNA 0.70–0.88. Protein engages its target.
- **Active site ordered:** catalytic-residue pLDDT (cat_min) 77–89 in every design.
- No design collapsed. AF3 finds them all structurally sound Cas9:sgRNA:DNA complexes.

## R-loop iPTM is NOT usable here (metric artifact, not a design failure)
sgRNA:DNA iPTM = 0.01–0.06 for designs **and** the ncas9 calibrator. Two reasons:
1. the 98 nt sgRNA is mostly scaffold, so the whole-chain iPTM is diluted by the non-pairing region;
2. AF3 keeps the DNA partly **duplexed** (DNA:DNA pairs 0.24–0.33 > sgRNA:DNA 0.01–0.03), i.e. it
   does not form the R-loop. → Do not rank on R-loop iPTM; defer R-loop/PAM to wet lab.

## Usable secondary signal: Cas9:sgRNA iPTM + pLDDT (ranked)
| Tier | designs |
|---|---|
| STRONG (gRNA≥0.79, pLDDT≥84) | d00107, d00824, d00523, d00526, d00097, d00900*, d00941, d00806* |
| ok | d00018, d00918 |
| **weaker — deprioritise** | d00272 (pLDDT 80.7), d00542* (gRNA 0.75), d00134 (0.75), d00101 (gRNA 0.70, pLDDT 81.6) |

\*consensus designs. The 3 weak-gRNA designs (d00101, d00134 bold-novelty + d00542) are among the
most divergent — weaker sgRNA engagement at high divergence is plausible.

## Gaps / caveats
- **Calibration incomplete:** only bind_dead + ncas9 folded; the critical **WT (positive)** and
  **scram_rec (fold-negative)** were NOT folded — so the PASS band isn't anchored, and we can't yet
  confirm AF3 even penalises nonsense for this scaffold.
- **bind_dead didn't discriminate:** its bridge-helix Arg→Glu still shows high binding (Cas9:gRNA 0.88,
  Cas9:DNA 0.95) — AF3 doesn't register that disruption.

## Recommendation
1. **Fold the 3 missing calibrators** (wt_spry, dcas9, scram_rec) — esp. scram_rec (fold floor) +
   wt_spry (confirms R-loop is an artifact). 3 cheap jobs; completes the gate.
2. **Carry the STRONG + consensus designs forward**; keep 1–2 bold ones (d00097/d00107) for the
   divergence test; deprioritise d00101/d00272.
3. **Wet-lab PAM panel remains the real arbiter** — AF3 cannot test catalysis or PAM relaxation.
