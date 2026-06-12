"""MSA parsing / encoding shared by conservation and coevolution stages."""
from __future__ import annotations
import gzip
import numpy as np

# 20 AAs + gap. Order matters: keep gap last so AA_GAP == 20.
AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY-"
AA_TO_INT = {aa: i for i, aa in enumerate(AA_ALPHABET)}
N_ALPHA = len(AA_ALPHABET)          # 21
AA_GAP = AA_TO_INT["-"]             # 20


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def parse_msa(fasta_path):
    """Return list of (id, sequence) preserving file order (first = WT query)."""
    seqs, cur_id, cur = [], None, []
    with _open(fasta_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur_id is not None:
                    seqs.append((cur_id, "".join(cur)))
                cur_id, cur = line[1:].split()[0], []
            else:
                cur.append(line.strip())
        if cur_id is not None:
            seqs.append((cur_id, "".join(cur)))
    if not seqs:
        raise ValueError(f"No sequences parsed from {fasta_path}")
    L = len(seqs[0][1])
    for sid, s in seqs:
        if len(s) != L:
            raise ValueError(f"MSA not aligned: {sid} has length {len(s)} != {L}")
    return seqs


def encode_msa(seqs):
    """(N, L) uint8 array; unknown chars -> gap."""
    N, L = len(seqs), len(seqs[0][1])
    msa = np.full((N, L), AA_GAP, dtype=np.uint8)
    for i, (_, seq) in enumerate(seqs):
        for j, c in enumerate(seq):
            msa[i, j] = AA_TO_INT.get(c.upper(), AA_GAP)
    return msa


def filter_columns(msa, max_gap_frac=0.5):
    """Drop columns with > max_gap_frac gaps. Returns (filtered_msa, kept_col_idx)."""
    gap_frac = np.mean(msa == AA_GAP, axis=0)
    keep = gap_frac <= max_gap_frac
    return msa[:, keep], np.where(keep)[0]
