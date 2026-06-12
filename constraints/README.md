# `constraints/` — Cas9 design-time constraint stage

Turns a **family MSA + target structure** into the fixed-position lists fed to
LigandMPNN. Formalizes the previously ad-hoc `/tmp` scripts into a permanent,
tested pipeline stage (Task #41).

## Two complementary signals

| | Conservation | Coevolution (DCA) |
|---|---|---|
| Module | `conservation.py` | `coevolution.py` |
| Order | 1st (single column) | 2nd (column pairs) |
| Question | "is THIS position fixed?" | "are these TWO coupled?" |
| Output | `top_frac` per position | per-position coupling score |
| Role in tiers | backbone (T2 ≥90%, T3 ≥50%) | +N positions conservation misses |

On the audited Cas9 set: **35 of the top 42 DCA positions are invisible to
conservation ≥50%** — that is the unique signal this stage adds.

## DCA backends (`--method`)
- **`mi_apc`** (default) — mutual information + average-product correction. No
  external deps, deterministic. This is the validated method behind the audited
  42-position set.
- **`evcouplings`** — mean-field DCA (CN/APC score) via the `evcouplings`
  package, focus mode. More rigorous global model. ⚠️ Mean-field inverts an
  `(L·20)²` covariance matrix, so on full-length Cas9 (L≈1341 → ~27k×27k, ~6 GB)
  it is memory-heavy and slow; best for focused sub-regions or a high-RAM box.
  Use `mi_apc` (which scales fine) as the default; switch to `evcouplings` when
  you want the global model on a region. (plmDCA via `plmc` can be slotted in
  later; `plmc` is not currently installed.)

## The mapping bug this stage fixes
`mapping.py` maps **MSA column → PDB author residue number** gap-safely. The old
`/tmp` script numbered positions 1,2,3… by ungapped index, but 8SRS chain A
starts at author residue **3** with internal gaps — every DCA position was
mis-assigned. The fix: the i-th ungapped WT residue → the i-th *present* author
residue. Validated: **0 non-SpRY WT/structure mismatches** across all 1341 columns.

## Usage
```bash
python -m constraints.cli \
    --msa constraints/data/cas9_msa.fasta.gz \
    --structure inputs/8SRS.cif \
    --method mi_apc \
    --n-top-dca 42 \
    --emit-fixed \
    --out outputs/constraints/8srs_constraints.json
```
Outputs per-position conservation tiers, the top-N DCA positions (PDB-numbered),
the DCA-unique-vs-conservation count, and (with `--emit-fixed`) a ready
`--fixed_residues` string (conserved ≥50% ∪ DCA ∪ SpRY).

## Tests
```bash
python -m pytest tests/ -q                  # fast: mapping + conservation
RUN_SLOW_DCA=1 python -m pytest tests/ -q   # + full DCA reproduction of the 42
```

## Files
```
msa_io.py        MSA parse/encode/gap-filter (handles .gz)
mapping.py       gap-safe MSA-column ↔ PDB author-number mapping
conservation.py  1st-order per-column conservation
coevolution.py   2nd-order DCA: mi_apc + evcouplings backends
pipeline.py      orchestrator + tier-merge helpers
cli.py           command-line entry point
data/cas9_msa.fasta.gz   3,782-seq Cas9 family MSA (row 0 = SpCas9 WT query)
```
