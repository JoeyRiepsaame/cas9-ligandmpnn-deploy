# SpRY-Cas9 Path-4.1b — Pipeline Flow

**Goal:** reverse-engineer a patentable SpRY-Cas9 — diverge enough for novelty (<70% identity
where safe) while preserving function (cutting + the relaxed-PAM phenotype). Strategy is
**per-domain**: PROTECT the PI/PAM domain + HNH; REDESIGN the REC lobe; MODERATE the RuvC/BH scaffold.

Legend: ✅ done & validated · 🔨 partially built · ⏳ pending

```
┌─ PHASE 0 · INPUTS ───────────────────────────────────────────────────────────┐
│  inputs/8SRS.cif            constraints/data/        ref/Q99ZW2_SpCas9         │
│  SpRY-Cas9 complex          cas9_msa.fasta.gz        _canonical.fasta          │
│  chain A=Cas9(1341)         3,782 Cas9 orthologs     UniProt, 1368 aa          │
│  B=gRNA C/c/D=DNA + Mg²⁺    (WT query = row 0)       (numbering ground truth)  │
└──────────┬──────────────────────────┬──────────────────────┬──────────────────┘
           ▼                          ▼                       │
┌─ PHASE 1 · CONSTRAINT FOUNDATION ──────────────────────┐    │
│  constraints/ stage (Task #41, pre-existing)            │    │
│   conservation (1st-order)  +  DCA mi_apc (2nd-order)   │    │
│   → outputs/constraints/8srs_constraints.json   ✅      │    │
│       conservation_tiers ge_90/70/50 · dca_top(42)      │    │
│       · spry_positions(11)                              │    │
│  domain_map.py → domain_map.json                ✅      │    │
│   PI(1099-1368)+HNH(775-908) = PROTECT                  │    │
│   REC1/2/3(94-717) = REDESIGN (novelty)                 │    │
│   RuvC-I/II/III+BH = MODERATE                           │    │
│   anchors: 11 SpRY (10 in PI) · 7 catalytic             │    │
└───────────────────────────┬─────────────────────────────┘    │
                            ▼                                  │
┌─ PHASE 2 · TIER CONSTRUCTION ──────────────────────────┐     │
│  build_tiers.py → tiers/*.txt + tiers_summary.json  ✅  │     │
│  always-fixed = SpRY(11)∪catalytic(7)∪contacts(35)      │     │
│                 ∪DCA(42) = 95   + per-domain cons floor: │     │
│   ┌────────────┬──────────┬───────────┬──────────┐      │     │
│   │            │aggressive│ BALANCED  │  safe    │      │     │
│   │ REC (redes)│  6-13%   │  23-50%   │ 41-67%   │      │     │
│   │ PI (protect)│   52%   │   79%     │  79%     │      │     │
│   │ HNH(protect)│   62%   │   84%     │  84%     │      │     │
│   │ TOTAL fixed │  31.4%  │  57.5% ◄  │ 69.8%    │      │     │
│   └────────────┴──────────┴───────────┴──────────┘      │     │
└──────────────┬──────────────────────────┬───────────────┘     │
               ▼                          ▼                      ▼
┌─ extract_wt.py ─────────┐   ┌─ PHASE 3 · VALIDATION GATES (pre-GPU) ──────────┐
│ spry_cas9_wt.fasta  ✅  │──►│ audit_numbering.py            ✅ PASS           │
│ wt_resnum_to_index.json │   │  A external: author# == canonical# (0 offsets)  │
│ (catalytic+SpRY audited)│   │  B internal: MSA-vs-struct diffs == 11 SpRY     │
└─────────────────────────┘   │  C pipeline: map 0-mismatch; cons/DCA mapped    │
                              │ pre_audit.py (blind-spot)     ✅ PASS           │
                              │  omit-C: only C80 forced (<50% cons); C574      │
                              │   protected → 0 correlation artifacts           │
                              │  bias: 0 (E-penalty dropped) → no tier change   │
                              └───────────────────────┬─────────────────────────┘
                                                      │
═══════════════════════ GPU BOUNDARY (first compute spend) ═══════════════════════
                                                      ▼
┌─ PHASE 4 · GENERATION ⏳ ───────────────────────────────────────────────────────┐
│  deploy_path4_1b.sh (preflight ✅; needs RunPod/MPS)                              │
│   per tier: LigandMPNN(atom ctx)+ProteinMPNN+SolubleMPNN                          │
│   --omit_AA C   ✗ NO --bias_AA   catalytic 7 = fixed WT                           │
│   ~510 designs/tier = ~1,530 total                                               │
│   → outputs/<tier>/.../seqs/*.fa + ranked_designs.csv (MPNN score = 1 axis)       │
└───────────────────────────────────────┬─────────────────────────────────────────┘
                                        ▼
┌─ PHASE 5 · SCORING & SELECTION ⏳ ──────────────────────────────────────────────┐
│  esm_score_windowed.py 🔨  windowed ESM-2 (1341>1022) naturalness                │
│  dutton_correct.py     ⏳  de-bias MPNN log-odds (principled E-correction)        │
│  pareto_rank.py        ⏳  Pareto frontier: ESM ⟂ Dutton-MPNN                     │
│  audit.py              ⏳  byte-match sequence-safety audit                       │
│   → per-tier frontier → operating tier → synthesis_shortlist.json                │
└───────────────────────────────────────┬─────────────────────────────────────────┘
                                        ▼
┌─ PHASE 6 · AF3 STRUCTURAL TRIAGE ⏳ (sanity gate, NOT a ranker) ─────────────────┐
│  build_af3 (8SRS layout) → batch · manual AlphaFold-Server upload                 │
│  parse_af3.py: chain0=Cas9, DNA frags, sgRNA; R-loop iPTM = cp[2][4]              │
│   R-loop iPTM + per-residue pLDDT (incl. ex-Cys site)                             │
└───────────────────────────────────────┬─────────────────────────────────────────┘
                                        ▼
┌─ PHASE 7 · THREE CATALYTIC TRACKS (post-hoc point mutations, ~0 extra compute) ──┐
│  active(7 cat WT) ─┐                                                              │
│  nickase H840A     ├─► point-mutate SELECTED winners (no re-generation)           │
│  dCas9 D10A+H840A ─┘                                                              │
└───────────────────────────────────────┬─────────────────────────────────────────┘
                                        ▼
┌─ PHASE 8 · WET-LAB PILOT ⏳ ────────────────────────────────────────────────────┐
│  the only true functional test (AF3 can't test PAM relaxation):                   │
│  cleavage on NGG + NAA/NAC/NAG PAMs → calibrate predictor-vs-function             │
│  before the large synthesis order                                                 │
└──────────────────────────────────────────────────────────────────────────────────┘

ALWAYS-FIXED ANCHORS carried through every stage: 11 SpRY · 7 catalytic
```

## How phylogeny enters the pipeline (verified against code)

Phylogeny is handled as a **confound-correction in the DCA stage only**, not as a positive design signal:

- **Henikoff reweighting** (`coevolution.py: sequence_weights(threshold=0.8)`) downweights
  redundant clades before coupling calc — **neff = 114 effective / 3,782 raw (~33× redundancy)**.
- **APC** removes residual phylogenetic/entropic background from the MI matrix.
- (`evcouplings` backend: `theta=0.2`, i.e. reweight at 80% identity — same principle.)

**Gaps (deliberately noted):**
- **Conservation tiers are phylogeny-NAIVE** — `conservation.py` uses raw unweighted `top_frac`.
  Since conservation is the fixation backbone, this is an internal inconsistency vs the DCA stage.
  *Recommended cheap fix:* apply the same Henikoff weights inside `compute_conservation`.
- **No RT-style explicit phylogenetics** ported to Cas9: no phylogenetic depth cutoff on the MSA,
  no pySCA/SCA sectors, no ancestral-sequence reconstruction, no comparative-genomics
  ancestral-state validation. (In the RT campaign these caught real artifacts, e.g. C90–184.)

## File index (`path4_1b/`)
| File | Role | Status |
|---|---|---|
| `domain_map.py` / `.json` | per-domain strategy + budget | ✅ |
| `build_tiers.py` / `tiers/` / `tiers_summary.json` | per-domain graded tiers | ✅ |
| `extract_wt.py` / `spry_cas9_wt.fasta` / `wt_resnum_to_index.json` | WT reference | ✅ |
| `audit_numbering.py` / `ref/Q99ZW2_*.fasta` | 3-way numbering proof | ✅ |
| `mpnn_blind_spot_detector.py` / `pre_audit.py` / `pre_audit_results.json` | pre-design artifact audit | ✅ |
| `esm_score_windowed.py` | windowed ESM-2 (the Cas9 length tweak) | 🔨 |
| `deploy_path4_1b.sh` | generation (no E-penalty, mixed MPNN) | 🔨 preflight ✅ |
| `dutton_correct.py` / `pareto_rank.py` / `audit.py` / AF3 builder+parser | post-generation | ⏳ |
```
```
