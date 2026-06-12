# SpRY-Cas9 (8SRS) LigandMPNN — Constraint Evidence Audit + Conservation vs DCA

**Date:** 2026-06-12
**Project:** SpRY-Cas9 3-tier LigandMPNN deploy (RunPod) — evidence-based constraint selection
**Session:** continues e79afab3-674a-42b3-a263-289983483cce
**Status:** Evidence analysis COMPLETE (8/8 ralph-loop iterations, 12 papers) → results audited & validated → deploy script pending

---

## 1. Problem / Goal

Before generating the 3-tier LigandMPNN deploy script for SpRY-Cas9 (8SRS), validate the evidence-based
constraint conclusions with an independent sanity check/audit, and clarify whether the current pipeline
implements **MSA Conservation** and how that differs from **DCA / coevolution**.

---

## 2. Audit Verdict — VALIDATED (with 2 flags)

Independent recount from source JSON (not trusting the narrative scripts). All headline numbers reconcile.

| Check | Result | Status |
|---|---|---|
| Tier counts vs /1341 residues | T1=53 (4.0%), T2=692 (51.6%), T3=1120 (83.5%) | ✅ exact |
| Tier nesting | T1 ⊂ T2 ⊂ T3 (no orphan positions) | ✅ clean |
| DCA top set | 42 positions, MI+APC cutoff 31.77 | ✅ |
| **DCA orthogonality** | **35 / 42 DCA positions NOT captured by conservation ≥50%** | ✅ confirms DCA adds unique signal |
| Same MSA for both signals | conservation + DCA both read `cas9_msa.fasta` (~3,782 seqs) | ✅ apples-to-apples |
| Literature aggregate 1/85 (1.2%) | = Sumida site 0/24 + King site 1/20 + King 10Å 0/32 + RamR 0/9 | ✅ arithmetic holds |
| All 6 cross-paper sanity checks | internally consistent across 7 studies | ✅ pass |

### Conservation tier sizes (PDB-numbered, /1341)
| Threshold | Positions | % of protein |
|---|---|---|
| ≥90% | 490 | 36.5% |
| ≥70% | 741 | 55.3% |
| ≥50% | 1048 | 78.2% |
| ≥25% | 1334 | **99.5%** ⚠️ |

### ⚠️ Audit Flag 1 — the 25% conservation cutoff is meaningless for OUR MSA
Tao et al. used 25% / 50% cutoffs on their **prime-editor RT** MSAs (more divergent). On our **deep Cas9
ortholog alignment (3,782 seqs)** a 25% threshold fixes 99.5% of positions — it is not a useful knob.
**Our meaningful conservation cutoffs are 50 / 70 / 90%.** Do NOT transplant Tao's 25% number directly;
the *direction* of their finding transfers, the *absolute threshold* does not.

### ⚠️ Audit Flag 2 — RamR is a different failure mode
The 1/85 "distance/site-only" denominator lumps in Clark-ElSayed's RamR (0/9), which fails for an
**allosteric/conformational-freezing** reason, not a distance-cutoff reason. Excluding it: **1/76 = 1.3%**
— conclusion unchanged (distance-only ≈ 1% vs conservation-based ≈ 29%, still a ~25–40× gap).

### Conclusion stands
- **Distance cutoffs are NOT critical** and can be counterproductive (King 10Å sphere = 0/32 active).
- **Conservation is the #1 determinant** (every paper).
- **DCA/coevolution is a real orthogonal signal** — empirically, 35/42 of our top DCA positions are invisible to conservation≥50%.
- **MPNN log-likelihood** = best post-hoc activity predictor (Johnson 2025) → add `--save_score 1`.
- Optimal fixed fraction sweet spot **54–73%** → our **T2 (52%)** and **T3 (84%)** bracket it; **T3 ≈ Tao PE8 winners (64–84%)**.

---

## 3. Does our pipeline implement MSA Conservation? How does it differ from DCA?

### Short answer
- **MSA Conservation: YES, implemented** — it's the 8SRS conservation analyzer output (`cas9_conservation.json`, per-column `top_frac` + `entropy` over the 3,782-seq MSA). This is what builds the conservation tiers.
- **DCA/Coevolution: ported but NOT yet a permanent pipeline stage** — it lives as the standalone `cas9_dca.py` (MI+APC), adapted from the RT/CmR projects. Formalizing it is the still-open **Task #41 ("add DCA/EVcouplings to Cas9 pipeline")**.

### The conceptual difference (this is the key point)

| | **MSA Conservation** | **DCA / Coevolution** |
|---|---|---|
| Order of statistic | **1st order** (single column) | **2nd order** (column *pairs*) |
| Question it answers | "Is *this* position fixed across evolution?" | "Do *these two* positions co-vary / are they coupled?" |
| Computation | per-column AA frequency: `top_frac` = freq of most common residue; `entropy` | mutual information between every column pair, **+ APC** (average-product correction) to remove phylogenetic/entropic background |
| What it captures | residues evolution won't touch (active site, fold core) | **epistatic networks** — pairs that must mutate *together*; allosteric/contact couplings |
| Blind spot | misses a residue that is variable *alone* but constrained *jointly* | needs a deep, diverse MSA; noisier per-position |
| In our data | ≥90%: 490 pos; ≥50%: 1048 pos | 42 top positions, **35 of which conservation≥50% misses** |

**Why both matter for Cas9:** conservation locks the obviously-essential core; DCA rescues the 35 "quietly
coupled" positions that a single-site view would have freed for redesign — exactly the multi-mutation
epistasis that AiCE (Fei et al. 2025, Cell) showed is essential for combinatorial designs to stay functional.
They are computed from the **same MSA** but extract **complementary** signal (1st- vs 2nd-order).

> Note: cas9-discovery's `homology_analysis.py` does **BLAST homology for novel-Cas9 discovery** — that is a
> *different* use of homology and is NOT the design-time conservation/DCA discussed here.

---

## 4. Final Tier Definitions (audited) — for the deploy script

Source: `/tmp/cas9_final_tiers.json` (validated this session). Target: SpRY-Cas9 8SRS chain A, 1341 residues.

| Tier | Fixed positions | % fixed | Composition | Predicted outcome |
|---|---|---|---|---|
| **T1** | 53 | 4.0% | active-site contacts + SpRY only | Expresses, **won't cut** — expression control (95% conf. it fails for activity) |
| **T2** | 692 | 51.6% | consv≥90% (490) + DCA (35 unique) + contacts + PTM + SpRY | Most novel; tests DCA hypothesis (50–70% conf. functional) |
| **T3** | 1120 | 83.5% | consv≥50% (1048) + DCA + contacts + PTM + SpRY | Safest; = Tao PE8 winner range (80–90% conf. functional) |

Nesting verified: **T1 ⊂ T2 ⊂ T3**.

---

## 5. Evidence Base (12 papers, wet-lab validated)

7 direct + 5 indirect. Full detail in memory `tao_mpnn_constraint_evidence.md`. Headline:
**Distance/site-only ≈ 1.2% active (1/85) vs conservation-based ≈ 29% (~100/350+).**

1. Sumida 2024 (JACS) — TEV: site-only 0/24; 50% conserved = best, 26× kcat/Km.
2. King 2024 (Angew) — tP4H: **10Å sphere = 0/32 (direct distance failure)**; 70% conserved = champion; GriE 27/32.
3. Tao 2026 (Nat Biotech) — PE RTs: conservation strong, distance 15/18/20Å no difference; PE8 winners 64–84% fixed.
4. AiCE / Fei 2025 (Cell) — 8 proteins, NO distance cutoffs; MPNN freq β≥0.8 + coevolution (LD≥0.5 + SCA≥90th). Coevolution critical for combos.
5. Ramírez-Sarmiento 2026 (bioRxiv) — PHL7 PETase: 6Å+conservation, 2/31 active (protein-specific sensitivity).
6. Clark-ElSayed 2025 (bioRxiv) — RamR: LigandMPNN 0/9 (allosteric freeze — MPNN failure mode; flag for Cas9 HNH).
7. Johnson 2025 (Nat Biotech) — 500+ designs: MPNN log-likelihood = best single activity predictor.
+ ChatGPT set: OpenCRISPR (131/209 Cas9s edit, OC-1 = 403 muts, catalytic core invariant), PAMmla, ProMEP, ancestral Cas9, NovaIscB/OMEGA.

---

## 6. Recommended deploy-script settings (next step)

```
--omit_AA "C"                      # avoid spurious cysteines
--bias_AA "E:-1.0"                 # mild glutamate down-bias
--ligand_mpnn_use_atom_context 1   # nucleic-acid/ligand atom context
--save_score 1 --save_probs 1      # MPNN log-likelihood ranking (Johnson 2025)
--fixed_residues "<tier positions>"
# 3 models × 2 temps (0.1/0.3) × 2 seeds × 10 seqs = 120 designs / tier
```
Rank designs by MPNN score within each tier before AF3 triage.

---

## 7. Next Steps
1. **Generate 3-tier RunPod LigandMPNN deploy script** (primary deliverable; tiers above).
2. **Task #41** — formalize DCA/EVcouplings as a permanent Cas9 pipeline stage (currently `cas9_dca.py` ad-hoc).
3. Consider EVcouplings as a more rigorous DCA replacement (plmDCA vs our MI+APC) when formalizing.
4. Post-design: MPNN-score rank → AF3 → order synthesis panel.

## 8. Key Files
- `/tmp/cas9_final_tiers.json` — audited tier definitions
- `/tmp/cas9_conservation.json` — per-column conservation (1st-order)
- `/tmp/cas9_dca_results_corrected.json` — 42 DCA positions, corrected PDB mapping (2nd-order)
- `/tmp/cas9_dca.py` — MI+APC DCA implementation (to be formalized, Task #41)
- `/tmp/audit_numbers.json` — this session's recount block
- Memory: `tao_mpnn_constraint_evidence.md`
