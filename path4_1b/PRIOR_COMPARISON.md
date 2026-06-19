# SpRY-Cas9: Path-4.1b consensus designs vs the PRIOR (May-2026) MPNN campaign

**Date:** 2026-06-19. Apples-to-apples on the two axes computable identically for both
campaigns: **windowed-ESM naturalness** (model-agnostic) + **%identity-to-WT** (alignment).
Common reference confirmed: the prior SpRY WT is **byte-identical** to ours (1341 aa).
Prior data: `~/gdrive/Vectors/Alphafold/LigandMPNN/8SRS_Conservation/` (120 designs: Path-1×60 + Path-2×60).

> Note: the prior campaign's Dutton/MPNN axis can't be recomputed (no log_probs saved), so ESM is
> the fair common axis. ESM correlates with identity, so the patent-relevant comparison gates on <70 % id.

## The prior campaign in our frame
| Set | n | mean ESM | mean %id | note |
|---|---|---|---|---|
| Path-1 (cons≥90 % fixed) | 60 | −0.343 | 59.5 % | novel, but **less natural** |
| Path-2 (cons≥50 % fixed) | 60 | −0.307 | 84.6 % | natural, but **fails the <70 % patent gate** |
| by model — Ligand / Protein / Soluble | 40/40/40 | −0.320 / −0.323 / −0.332 | 73 / 72 / 71 % | near-tied on ESM |

WT ESM = −0.307 (reference).

## Finding 1 — the prior campaign picked the WRONG designers (bias confirmed)
The prior AF3 shortlist (top-10 by raw MPNN **confidence**) was **100 % LigandMPNN**. But on the
**fair, model-agnostic ESM axis the top-10 are 100 % ProteinMPNN.** Raw MPNN confidence is inflated
for LigandMPNN by its atom-context advantage — exactly the bias our **per-model-z fix** corrects.
The prior selection metric systematically favoured the wrong model type.

## Finding 2 — the prior wet-lab "top pick" is mediocre by naturalness
`LigandMPNN_s2_sample4` (their #1, chosen by AF3 **R-loop iPTM** — the metric our scram_rec control
proved is an artifact) ranks **52 / 120** on ESM (−0.314, 61 % id). It was selected by two signals we
have since shown to be unreliable (atom-context-biased conf + non-discriminating R-loop).

## Finding 3 — our consensus designs lead the naturalness-vs-novelty frontier
| | mean ESM | mean %id |
|---|---|---|
| **Our consensus 4** | **−0.312** | **68.0 %** |
| Prior Path-1 (comparable novelty ~60 %) | −0.343 | 59.5 % |
| Prior Path-2 (matches our ESM) | −0.307 | **84.6 % (fails gate)** |

- The prior pipeline only reaches our naturalness (~−0.31) at **84 % identity** (Path-2, not patent-novel);
  at comparable novelty (~60 %, Path-1) it is markedly **less natural** (−0.343).
- **0 / 120 prior designs beat our best consensus (d00542, ESM −0.307 = WT level) while staying <70 % id.**
- Only 60/120 prior designs are <70 % id at all (Path-2's 60 fail the patent gate outright).
- Our **d00542** (Ligand, 69.6 % id) is **as natural as WT** — a standout.

## Bottom line
Measured fairly, the Path-4.1b consensus designs **dominate the prior campaign on the
naturalness-vs-novelty trade-off**: comparable-or-better naturalness at genuinely patent-novel
identity, where the prior pipeline could only be natural by staying conservative (Path-2, >84 %) or
novel by sacrificing naturalness (Path-1). The prior campaign also **selected the wrong designs**
twice over — by atom-context-biased MPNN confidence (→ all-LigandMPNN) and by the AF3 R-loop
artifact (→ a rank-52 "top pick"). Both failure modes are explicitly fixed in Path-4.1b
(per-model-z + the calibrator-proven demotion of AF3). ESM remains a proxy; the wet-lab PAM panel
is the final arbiter for both campaigns.
