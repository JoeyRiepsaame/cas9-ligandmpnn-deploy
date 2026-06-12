"""
MSA conservation — 1st-order, single-column signal.

Answers: "is THIS position fixed across the family?" Per-column frequency of the
most common residue (gaps excluded), plus Shannon entropy. This is the
conservation backbone of the LigandMPNN tiers (T2 = consv>=90%, T3 = consv>=50%).
"""
from __future__ import annotations
import numpy as np
from .msa_io import parse_msa, encode_msa, AA_GAP, N_ALPHA, AA_TO_INT

_INT_TO_AA = {v: k for k, v in AA_TO_INT.items()}


def compute_conservation(fasta_path):
    """
    Returns list (len = alignment length) of per-column dicts indexed by ORIGINAL
    column: {col, query_aa, top_aa, top_frac, n_seqs, entropy}.
    top_frac excludes gaps from the denominator (consensus among aligned seqs).
    Row 0 is the WT query.
    """
    seqs = parse_msa(fasta_path)
    msa = encode_msa(seqs)
    _, L = msa.shape
    out = []
    for col in range(L):
        column = msa[:, col]
        non_gap = column[column != AA_GAP]
        n = int(non_gap.size)
        query_aa = _INT_TO_AA[int(msa[0, col])]
        if n == 0:
            out.append({"col": col, "query_aa": query_aa, "top_aa": "-",
                        "top_frac": 0.0, "n_seqs": 0, "entropy": 0.0})
            continue
        counts = np.bincount(non_gap, minlength=N_ALPHA).astype(float)
        freqs = counts / n
        top_i = int(np.argmax(counts))
        nz = freqs[freqs > 0]
        out.append({
            "col": col,
            "query_aa": query_aa,
            "top_aa": _INT_TO_AA[top_i],
            "top_frac": round(float(freqs[top_i]), 4),
            "n_seqs": n,
            "entropy": round(float(-(nz * np.log2(nz)).sum()), 4),
        })
    return out


def conserved_positions(conservation, col_to_pdb, threshold):
    """PDB author resnums whose (original-column) top_frac >= threshold."""
    keep = set()
    for col, pdb_num in col_to_pdb.items():
        if conservation[col]["top_frac"] >= threshold:
            keep.add(pdb_num)
    return keep
