"""
Constraint pipeline orchestrator.

Permanent Cas9 design stage: MSA + structure -> conservation + coevolution(DCA),
both mapped to PDB author numbering, merged into LigandMPNN fixed-position tiers.

  run_constraints()      end-to-end; writes a results JSON
  merge_dca_into_tiers() add the top-N coupled positions to an existing tier set
"""
from __future__ import annotations
import json
import numpy as np

from .msa_io import parse_msa
from .conservation import compute_conservation, conserved_positions
from .coevolution import compute_coevolution, top_coupled_positions
from .mapping import parse_structure_chain, build_column_to_pdb

# SpRY positions (SpCas9 WT query != SpRY structure) — mismatches here are expected,
# and these are always fixed (engineered PAM-relaxing substitutions must be kept).
SPRY_POSITIONS = [61, 1111, 1135, 1136, 1218, 1219, 1317, 1322, 1333, 1335, 1337]


def run_constraints(msa_path, structure_path, chain="A",
                    method="mi_apc", n_top_dca=42,
                    conservation_thresholds=(0.90, 0.70, 0.50),
                    fixed_for_dca_exclusion=0.90,
                    out_json=None, **dca_kwargs):
    """
    Compute conservation + coevolution and map both to PDB numbering.

    fixed_for_dca_exclusion: DCA's value is the positions conservation MISSES, so
        when picking top-N coupled positions we exclude those already conserved at
        this threshold (default 0.90 — i.e. report DCA hits beyond the strong core).

    Returns a results dict (also written to out_json if given).
    """
    seqs = parse_msa(msa_path)
    wt_aligned = seqs[0][1]

    structure = parse_structure_chain(structure_path, chain=chain)
    col_to_pdb, map_report = build_column_to_pdb(
        wt_aligned, structure, spry_positions=SPRY_POSITIONS, strict=True)

    conservation = compute_conservation(msa_path)

    cons_tiers = {f"ge_{int(t*100)}": sorted(conserved_positions(conservation, col_to_pdb, t))
                  for t in conservation_thresholds}

    coevo = compute_coevolution(msa_path, method=method, **dca_kwargs)
    strong_core = conserved_positions(conservation, col_to_pdb, fixed_for_dca_exclusion)
    dca = top_coupled_positions(coevo, col_to_pdb, fixed_already=strong_core,
                                n_top=n_top_dca, exclude_fixed=True)

    # how much unique signal does DCA add beyond conservation>=50%?
    cons50 = set(conserved_positions(conservation, col_to_pdb, 0.50))
    dca_unique_vs_50 = sorted(set(dca["positions"]) - cons50)

    results = {
        "inputs": {"msa": str(msa_path), "structure": str(structure_path),
                   "chain": chain, "method": method},
        "mapping_report": {k: v for k, v in map_report.items()
                           if k != "unexpected_mismatches"} |
                          {"n_unexpected_mismatches": len(map_report["unexpected_mismatches"])},
        "coevolution_meta": {k: coevo[k] for k in
                             ("method", "n_sequences", "alignment_length",
                              "filtered_columns", "neff")},
        "conservation_tier_sizes": {k: len(v) for k, v in cons_tiers.items()},
        "conservation_tiers": cons_tiers,
        "dca_top": dca["positions"],
        "dca_score_cutoff": dca["score_cutoff"],
        "dca_ranked": dca["ranked"],
        "dca_unique_vs_consv50": dca_unique_vs_50,
        "n_dca_unique_vs_consv50": len(dca_unique_vs_50),
        "spry_positions": SPRY_POSITIONS,
    }
    if out_json:
        with open(out_json, "w") as f:
            json.dump(results, f, indent=2)
    return results


def merge_dca_into_tiers(tier_positions, dca_positions):
    """Union an existing tier fixed-list with DCA positions; return sorted ints."""
    return sorted(set(int(x) for x in tier_positions) | set(int(x) for x in dca_positions))


def fixed_residues_string(positions, chain="A"):
    """LigandMPNN --fixed_residues token string, e.g. 'A53 A61 ...'."""
    return " ".join(f"{chain}{int(p)}" for p in sorted(set(positions)))
