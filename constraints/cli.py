#!/usr/bin/env python3
"""
CLI for the Cas9 constraint stage.

    python -m constraints.cli \
        --msa constraints/data/cas9_msa.fasta.gz \
        --structure inputs/8SRS.cif \
        --method mi_apc \
        --n-top-dca 42 \
        --out outputs/constraints/8srs_constraints.json

Backends: mi_apc (default, no deps) | evcouplings (mean-field DCA, needs evcouplings).
"""
from __future__ import annotations
import argparse
import json
import os
from .pipeline import run_constraints, fixed_residues_string


def main(argv=None):
    ap = argparse.ArgumentParser(description="Cas9 conservation + DCA constraint stage")
    ap.add_argument("--msa", required=True, help="family MSA (.fasta or .fasta.gz); row 0 = WT query")
    ap.add_argument("--structure", required=True, help="target structure (.cif preferred, or .pdb)")
    ap.add_argument("--chain", default="A")
    ap.add_argument("--method", default="mi_apc", choices=["mi_apc", "evcouplings"])
    ap.add_argument("--n-top-dca", type=int, default=42)
    ap.add_argument("--dca-exclusion-threshold", type=float, default=0.90,
                    help="exclude positions already conserved at this frac when ranking DCA")
    ap.add_argument("--out", default="outputs/constraints/constraints.json")
    ap.add_argument("--emit-fixed", action="store_true",
                    help="also write a fixed_residues string of conserved>=50%% + DCA")
    args = ap.parse_args(argv)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    res = run_constraints(
        args.msa, args.structure, chain=args.chain, method=args.method,
        n_top_dca=args.n_top_dca, fixed_for_dca_exclusion=args.dca_exclusion_threshold,
        out_json=args.out,
    )

    print("=== Cas9 constraint stage ===")
    print(f"method            : {res['inputs']['method']}")
    print(f"sequences         : {res['coevolution_meta']['n_sequences']}")
    print(f"filtered columns  : {res['coevolution_meta']['filtered_columns']}")
    print(f"mapping mismatches: {res['mapping_report']['n_unexpected_mismatches']} (expect 0)")
    print("conservation tiers:")
    for k, v in res["conservation_tier_sizes"].items():
        print(f"   {k:8} {v}")
    print(f"DCA top-{args.n_top_dca}        : {len(res['dca_top'])} positions, "
          f"cutoff {res['dca_score_cutoff']:.3f}")
    print(f"DCA unique vs c50 : {res['n_dca_unique_vs_consv50']} "
          f"(positions conservation>=50% misses)")
    print(f"written           : {args.out}")

    if args.emit_fixed:
        cons50 = set(res["conservation_tiers"]["ge_50"])
        merged = sorted(cons50 | set(res["dca_top"]) | set(res["spry_positions"]))
        s = fixed_residues_string(merged, chain=args.chain)
        fp = os.path.splitext(args.out)[0] + "_fixed.txt"
        with open(fp, "w") as f:
            f.write(s + "\n")
        print(f"fixed_residues    : {len(merged)} positions -> {fp}")


if __name__ == "__main__":
    main()
