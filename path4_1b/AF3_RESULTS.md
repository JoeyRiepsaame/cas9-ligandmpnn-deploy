# SpRY-Cas9 Path-4.1b — AF3 Triage Results (FINAL, 22 folds)

**Date:** 2026-06-19 · **18 designs + 4 calibrators** (bind_dead, dcas9, ncas9, **scram_rec**).
Table: `af3_interface_table.csv`. The full calibration ladder is now complete.

## HEADLINE: AF3 is NON-DISCRIMINATING for this campaign (the calibrator ladder proved it)
The fold-negative control **`scram_rec`** (a 60-residue stretch of REC2 randomly **shuffled** —
deliberate nonsense) scores **as high as the real designs globally**: iPTM 0.92, pLDDT 86.0,
Cas9:sgRNA 0.87, cat_min 82.8. Local check in the shuffled window (author 200–260):

| job | pLDDT @200–260 | global pLDDT |
|---|---|---|
| scram_rec (nonsense) | **43.7** | 89.6 |
| dCas9 (native REC2) | **88.6** | 94.2 |
| d00806 / d00900 / d00781 / d00107 (real designs) | 41.8 / 40.7 / 47.5 / 35.6 | ~84–90 |

→ The real redesigns score **as low as the scramble** in the redesigned REC2; only the **native**
sequence folds it confidently. **AF3 cannot distinguish a genuine heavy REC redesign from random
nonsense.** High global scores are carried by the conserved ~95 % of the protein. R-loop iPTM is
also an artifact (dCas9 binds perfectly yet R-loop = 0.03; 98 nt sgRNA dilution + AF3 keeps DNA duplexed).

## What AF3 DID establish (a weak sanity gate only)
- No design **globally collapses**: iPTM 0.90–0.92, Cas9:DNA 0.93–0.95, catalytic site ordered.
- That is the limit of its value here — it cannot **rank** designs or **validate the REC redesign**.

## Implication for selection
- **Do NOT use AF3 scores to narrow the pilot** — they don't discriminate (scram_rec passes).
- Selection rests on the **sequence-based axes** (ESM naturalness + per-model-z Dutton) already used,
  and the **wet-lab PAM panel is now the SOLE functional arbiter**.
- This is the RamR/AF3-blindness failure mode, empirically confirmed for Cas9 — exactly what the
  calibration ladder was built to detect. (Had we trusted "all designs pass," we'd have been misled.)

## Carry-forward (by sequence axes, AF3 neither confirms nor refutes)
- **Consensus designs** (in both shortlists) remain the safest: d00542, **d00781, d00806, d00900**.
- Span the wet-lab pilot across provenance/model/tier as planned (`pilot_set.json`).
- Wet-lab readout: cleavage on NGG + NAA/NAC/NAG.
