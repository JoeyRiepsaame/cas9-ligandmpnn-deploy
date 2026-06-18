# SpRY-Cas9 Path-4.1b — Wet-lab Pilot Plan (Phase 8)

Why a pilot: the <70%-identity frontier for a 1341-aa multidomain nuclease is
**unproven** (lit review: RT ~70–75%, single-domain down to ~41%). AF3 is a
fold/binding gate only — it cannot test catalysis or PAM relaxation. The pilot
calibrates predictor-vs-function before any large synthesis order.

## Constructs (pilot library)
From `synthesis_shortlist.json` (17 audited, <70% id), after AF3 triage:
- **~6–8 designs spanning both flavors**: T_aggressive (~57–61% id, bold) + T_balanced (~69%, safer),
  prioritising AF3 R-loop in the WT/dCas9 band + catalytic/ex-Cys pLDDT ≥70.
- **Controls (essential):**
  - WT SpRY-Cas9 (positive, relaxed-PAM active)
  - dCas9 (D10A+H840A) — no-cut control
  - a near-WT high-identity design (cleavage-likely sanity anchor)
- Codon-optimise for the expression host; **run the file-based sequence audit before ordering**.

## Assay
1. Express + purify (or IVTT) each Cas9 variant; reconstitute RNP with the 8SRS-matched sgRNA.
2. **In-vitro cleavage** on a dsDNA substrate panel that varies ONLY the PAM:
   - **NGG** (canonical, positive baseline) · **NAA · NAC · NAG** (the relaxed-PAM phenotype SpRY enables)
   - optional: NGA, NGC, NCG to map the PAM-relaxation profile.
3. Readout: % cleavage (gel or capillary) per PAM, vs WT SpRY and dCas9.

## What success looks like
- A design **cuts on NGG AND on ≥1 non-NGG PAM** at a meaningful fraction of WT SpRY
  → relaxed-PAM phenotype preserved at <70% identity (the patentable result).
- Designs that cut NGG only (lost PAM relaxation) → PI-domain redesign went too far;
  fall back to higher PI fixation.
- Designs dead on all PAMs → fold/active-site lost; correlate with AF3 R-loop/pLDDT to
  refine the predictor.

## Calibration / iteration
- Plot measured activity vs ESM, vs Dutton-MPNN, vs %identity, vs AF3 R-loop → learn which
  predictor (if any) tracks function at this divergence; **validate the predictor on this pilot
  before scaling** (ftMLDE warning — zero-shot predictors mis-rank at high divergence).
- Best PAM-relaxed cutters → derive nickase (H840A) / dCas9 (D10A+H840A) by post-hoc point
  mutation for base/prime-editing / CRISPRi applications.
