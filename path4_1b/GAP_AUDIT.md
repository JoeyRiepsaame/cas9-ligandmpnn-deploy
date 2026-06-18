# RT → SpRY-Cas9 Pipeline Gap Audit

Evidence-based comparison of the MMLV-RT Path-4.1b pipeline (the validated benchmark)
against the SpRY-Cas9 port, by inventorying **both** codebases (not memory).
RT source: `rt_sweep_20260501_154722/{scripts,sweep_results_4_1b,af3_inputs}`.
Date: 2026-06-18.

Buckets: ✅ ported & validated · ⏳ planned (post-generation, tracked) · ❗ missing / not yet planned

---

## A. Ported & validated ✅
| RT feature | Cas9 status |
|---|---|
| Graded fixation tiers | ✅ as **per-domain** tiers (`build_tiers.py`) |
| `--omit_AA C` | ✅ (`deploy_path4_1b.sh`); pre-audited (only C80 forced) |
| Drop E-penalty (Dutton replaces it) | ✅ no `--bias_AA` in deploy |
| ESM naturalness (COMPSS) | ✅ + Cas9 length tweak: **windowed** ESM-2 (`esm_score_windowed.py`) |
| Blind-spot detector + pre-design scan | ✅ (`mpnn_blind_spot_detector.py`, `pre_audit.py`) |
| Henikoff reweighting in **DCA** | ✅ already in `constraints/coevolution.py` (neff 114/3782) |
| Henikoff reweighting in **conservation** | ✅ **NEW — added in this port** (`reweight_conservation.py`); RT relied on pySCA instead |
| Numbering correctness | ✅ **stronger than RT** — 3-way audit vs canonical SpCas9 (`audit_numbering.py`) |

## B. Planned — pending generation output ⏳ (tracked, not missing)
These are faithful ports that simply need MPNN designs to run against:
| RT script | Cas9 plan |
|---|---|
| `extract_seqs.py` (dedup → unique_meta, n_copies) | ⏳ port |
| MPNN single-AA log-odds (`8WUV.pt`) | ⏳ emit `8SRS.pt` (needed for Dutton) |
| `dutton_correct.py` (de-bias MPNN log-odds) | ⏳ port |
| `pareto_rank.py` (ESM ⟂ MPNN + %id via Bio.Align) | ⏳ port |
| `audit.py` (reconciliation + byte-match seq-safety) | ⏳ port |
| `build_af3.py` / `parse_af3.py` | ⏳ adapt to **8SRS layout** (chain0=Cas9, DNA frags, sgRNA; R-loop=cp[2][4]) |
| `analyze.py` / `master_table` / `per_tier_rollup` | ⏳ reporting port |

## C. Missing / not yet planned ❗ (the substantive findings)

| # | RT feature | Why it mattered in RT | Cas9 gap | Priority |
|---|---|---|---|---|
| C1 | Conservation reweighting | phylogeny bias in fixation backbone | ✅ **CLOSED** (`reweight_conservation.py`) | — |
| C2 | **Design-ensemble LD / coevolution** (`ld_p35_within.json`, `ld_scoring_table.json`) + LD-compatibility scoring | found which positions co-vary *across designs* (distinct from MSA-DCA); scored designs for internal compatibility; caught the C90–184 artifact | ✅ **CLOSED** — `ld_analysis.py` (vectorised MI+APC, self-tested; 10K×600 in 2.7s) + joint-freq scoring; runs post-generation | — |
| C3 | **10K-ensemble** single-constraint run for robust LD/MI | 60-design LD over-called couplings; 10K gave stable stats | ✅ **script built** (`deploy_10k_ensemble.sh`, preflighted) — needs GPU run | (GPU) |
| C4 | **AF3 reference-panel calibration** | known-active (PE8d) + known-dead designs calibrated whether AF3 discriminates | **no Cas9 active/dead calibrators chosen** | **HIGH** |
| C5 | pySCA evolutionary coupling (Pfam) ∪ ensemble LD | added independent evolutionary-coupling signal | not ported (Cas9 RuvC/HNH Pfam exist) | MEDIUM |
| C6 | Comparative-genomics ancestral validation (412K BV-BRC) | validated ancestral states / refuted artifacts (I90) | not ported | MEDIUM |
| C7 | Negative-design controls (unfix suspected load-bearers) | lit-review upgrade #5 | not implemented as a tier variant | MEDIUM |
| C8 | Phylogenetic depth-cutoff curation of the MSA | RT chose ≥28%-to-MMLV w/ retroviral landmarks | Cas9 MSA implicitly 47–100% (median 60%); no explicit curation decision | MEDIUM (partly characterized now) |
| C9 | blind-spot `--msa` evolutionary validation pass | confirmed which artifacts have natural support | ran pre-scan only (0 artifacts predicted → low value) | LOW |
| C10 | GRACE explicit cavity/pocket fixing | lit-review upgrade #4 | partially covered by the 35 contacts; no dedicated pocket analysis | LOW |
| C11 | `ligandmpnn_runner.py` persistent runner | gradio/native ensemble runner | functionally replaced by `deploy_path4_1b.sh` (native) | none (equivalent) |
| C12 | position-specific tracking report (`position_409_report`) | tracked the flexible C409→E choice | optional analog for a Cas9 position of interest | LOW |
| C13 | predictor-vs-function pilot validation (lit warning B) | calibrate zero-shot predictors before big synthesis | Phase-8 wet-lab plan (conceptual) | (wet-lab) |

---

## Headline + recommended actions

**The pipeline is methodologically faithful on the *per-design selection* path** (tiers, omit-C,
no-E-penalty, ESM, Dutton, Pareto, blind-spot) — bucket B is just deferred-until-generation porting.

**The real gap is the *coevolution/ensemble* path (C2+C3) and AF3 calibration (C4)** — and C3 is
explicitly inside the "full Path-4.1b" scope the user chose, but the current `deploy_path4_1b.sh`
under-provisions it (~510/tier, no dedicated 10K LD ensemble). Recommended before/with generation:

1. **Add a 10K single-constraint ensemble on T_balanced** (separate from the selection sweep) +
   port the design-ensemble LD analysis + LD-compatibility scoring (C2+C3). Reconciles Task #4's
   "10K-ensemble" title with the deploy.
2. **Choose Cas9 AF3 calibrators** (C4): WT SpRY (active) + a known catalytically-dead/!cut variant
   as the known-negative, so the AF3 R-loop iPTM gate is interpretable.
3. Defer C5–C8 (pySCA, comparative genomics, negative-design, depth-curation) unless we want RT-level
   evolutionary rigor — they're rigor multipliers, not blockers.
