"""
Regression tests for the Cas9 constraint stage.

Fast tests (mapping + conservation) run in seconds and are always on.
The full MI+APC DCA reproduction is slow (~minutes); it runs only when the
reference file /tmp/cas9_dca_results_corrected.json is present AND
RUN_SLOW_DCA=1 is set, so CI stays fast.

    python -m pytest tests/ -q                 # fast only
    RUN_SLOW_DCA=1 python -m pytest tests/ -q  # include DCA reproduction
"""
import os
import sys
import json
import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)

from constraints.msa_io import parse_msa
from constraints.conservation import compute_conservation, conserved_positions
from constraints.mapping import parse_structure_chain, build_column_to_pdb
from constraints.pipeline import SPRY_POSITIONS, fixed_residues_string

MSA = os.path.join(HERE, "constraints/data/cas9_msa.fasta.gz")
STRUCT = os.path.join(HERE, "inputs/8SRS.cif")


@pytest.fixture(scope="module")
def mapping():
    seqs = parse_msa(MSA)
    structure = parse_structure_chain(STRUCT, "A")
    col_to_pdb, report = build_column_to_pdb(
        seqs[0][1], structure, spry_positions=SPRY_POSITIONS, strict=True)
    return col_to_pdb, report, seqs


def test_msa_loaded():
    seqs = parse_msa(MSA)
    assert len(seqs) == 3782
    assert seqs[0][0] == "SpCas9_WT_query"


def test_mapping_gapsafe_no_mismatches(mapping):
    col_to_pdb, report, _ = mapping
    # 8SRS chain A: 1341 modelled residues, author numbering 3..1366
    assert report["structure_residues"] == 1341
    assert report["structure_range"] == (3, 1366)
    assert report["columns_mapped"] == 1341
    # every WT residue must match the structure except the 11 SpRY substitutions
    assert len(report["unexpected_mismatches"]) == 0


def test_conservation_tier_sizes(mapping):
    col_to_pdb, _, _ = mapping
    cons = compute_conservation(MSA)
    sizes = {t: len(conserved_positions(cons, col_to_pdb, t)) for t in (0.90, 0.70, 0.50)}
    # audited values (2026-06-12)
    assert sizes[0.90] == 490
    assert sizes[0.50] == 1048
    assert sizes[0.70] == 741


def test_fixed_residue_string_format():
    s = fixed_residues_string([61, 10, 10, 1135], chain="A")
    assert s == "A10 A61 A1135"


def test_evcouplings_backend_runs():
    """EVcouplings mean-field DCA backend runs and returns sane scores.
    Windowed (400 cols / 1500 seqs) so the mean-field matrix inversion is small."""
    pytest.importorskip("evcouplings")
    import warnings
    import numpy as np
    from constraints.coevolution import _evcouplings_di
    seqs = parse_msa(MSA)
    sub = [(sid, s[:400]) for sid, s in seqs[:1500]]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ocs, focus = _evcouplings_di(sub, theta=0.2)
    vals = np.array(list(ocs.values()))
    assert len(ocs) == len(focus) > 0
    assert np.all(np.isfinite(vals)) and np.all(vals >= 0)


@pytest.mark.skipif(
    not (os.environ.get("RUN_SLOW_DCA") and
         os.path.exists("/tmp/cas9_dca_results_corrected.json")),
    reason="slow DCA reproduction; set RUN_SLOW_DCA=1 with reference file present",
)
def test_dca_reproduces_validated_42(mapping):
    from constraints.coevolution import compute_coevolution, top_coupled_positions
    col_to_pdb, _, _ = mapping
    saved = json.load(open("/tmp/cas9_dca_results_corrected.json"))
    path1b = set(int(x) for x in saved["path1b_fixed"])
    saved_top42 = sorted(int(x) for x in saved["top42_dca_positions"])

    coevo = compute_coevolution(MSA, method="mi_apc")
    dca = top_coupled_positions(coevo, col_to_pdb, fixed_already=path1b,
                                n_top=42, exclude_fixed=True)
    assert sorted(dca["positions"]) == saved_top42
