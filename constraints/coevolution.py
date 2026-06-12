"""
Coevolution / DCA — 2nd-order, column-pair signal.

Answers: "are these TWO positions coupled?" Captures the epistatic network that
single-site conservation is blind to. Formalizes the previously ad-hoc /tmp
cas9_dca.py (Task #41), with the position-mapping bug fixed (see mapping.py).

Two backends, selectable via `method`:
  "mi_apc"      : mutual information + average-product correction. No external
                  deps, fast, deterministic. This is the validated method that
                  produced the audited 42-position Cas9 set.
  "evcouplings" : mean-field DCA (DI score) via the evcouplings package. More
                  rigorous global model; used when available and requested.

Both return a per-column coupling score (sum of corrected couplings to all other
columns), which is then ranked to pick the top-N coupled, designable positions.
"""
from __future__ import annotations
import numpy as np
from .msa_io import parse_msa, encode_msa, filter_columns, N_ALPHA, AA_GAP


# ----------------------------------------------------------------------------- weights
def sequence_weights(msa, threshold=0.8, block=500):
    """Henikoff-style reweighting: w_n = 1/|seqs >= threshold identity to n|.
    Neff = sum(w). Deterministic. O(N^2 L) but vectorised per query."""
    N = msa.shape[0]
    w = np.ones(N, dtype=np.float64)
    for s in range(0, N, block):
        e = min(s + block, N)
        for i in range(s, e):
            ident = np.mean(msa == msa[i], axis=1)
            w[i] = 1.0 / max(int(np.sum(ident >= threshold)), 1)
    return w


# ----------------------------------------------------------------------------- MI+APC
def _mi_apc(msa, weights, pseudocount=0.5):
    """Per-column MI+APC coupling scores. Identical algorithm to the validated
    /tmp implementation, vectorised over amino-acid pairs."""
    N, L = msa.shape
    wsum = float(weights.sum())

    # single-site frequencies (with pseudocount), L x A.
    # NB: cast the boolean indicator to float64 BEFORE matmul — a float@bool
    # matmul under numpy 2.x can overflow/produce NaNs.
    f1 = np.zeros((L, N_ALPHA))
    for a in range(N_ALPHA):
        f1[:, a] = (weights @ (msa == a).astype(np.float64)) / wsum
    f1 = (1 - pseudocount) * f1 + pseudocount / N_ALPHA

    MI = np.zeros((L, L), dtype=np.float64)
    for a in range(N_ALPHA):
        Xa = (msa == a).astype(np.float64) * weights[:, None]    # N x L
        for b in range(N_ALPHA):
            Xb = (msa == b).astype(np.float64)                   # N x L
            fij = (Xa.T @ Xb) / wsum                             # L x L joint
            fij = (1 - pseudocount) * fij + pseudocount / (N_ALPHA * N_ALPHA)
            fmarg = np.outer(f1[:, a], f1[:, b])
            with np.errstate(divide="ignore", invalid="ignore"):
                contrib = fij * np.log(fij / (fmarg + 1e-10) + 1e-10)
            MI += np.nan_to_num(contrib, nan=0.0)
    np.fill_diagonal(MI, 0.0)

    # APC (average product correction)
    row = MI.mean(axis=1)
    apc = np.outer(row, row) / (MI.mean() + 1e-10)
    C = np.maximum(MI - apc, 0.0)
    np.fill_diagonal(C, 0.0)
    return C.sum(axis=1)


# ----------------------------------------------------------------------------- EVcouplings
def _evcouplings_di(seqs, theta=0.2):
    """
    Mean-field DCA via evcouplings, in FOCUS mode. Returns {orig_col: score}
    where score = sum of corrected (APC) DI to all other focus columns.

    Focus columns = the columns where the WT query (row 0) is non-gap; these are
    exactly the columns that map to a structure residue, so the result keys are
    directly the original MSA column indices. Raises ImportError if evcouplings
    is unavailable.
    """
    import io
    from evcouplings.couplings.mean_field import MeanFieldDCA
    from evcouplings.align import Alignment

    wt = seqs[0][1]
    focus_cols = [c for c, a in enumerate(wt) if a != "-"]          # WT non-gap columns
    n_focus = len(focus_cols)

    # Build a focus-mode a2m: row 0 header carries 1-based numbering so evcouplings
    # sets index_list = 1..n_focus. Focus columns uppercase; gaps as '-'.
    buf = io.StringIO()
    for k, (sid, seq) in enumerate(seqs):
        chars = []
        for c in focus_cols:
            ch = seq[c].upper()
            chars.append(ch if ch in "ACDEFGHIKLMNPQRSTVWY" else "-")
        header = f"query/1-{n_focus}" if k == 0 else f"{sid}/1-{n_focus}"
        buf.write(f">{header}\n{''.join(chars)}\n")
    buf.seek(0)

    aln = Alignment.from_file(buf, format="fasta")
    model = MeanFieldDCA(aln).fit(theta=theta)     # reweighting at identity 1-theta
    ecs = model.ecs                                # DataFrame: i, j, fn, cn (APC-corrected)

    # column coupling score = sum of positive CN (corrected-norm) couplings per
    # focus position. model.index_list is 1..n_focus; map back to original columns.
    score_by_focuspos = {}
    for _, row in ecs.iterrows():
        cn = float(row["cn"])
        if cn <= 0:
            continue
        for pos in (int(row["i"]), int(row["j"])):
            score_by_focuspos[pos] = score_by_focuspos.get(pos, 0.0) + cn

    # focus position p (1-based, in index_list) -> original column focus_cols[p-1]
    orig_col_score = {}
    for p, sc in score_by_focuspos.items():
        orig_col_score[int(focus_cols[p - 1])] = sc
    # ensure every focus column present (zero if no positive coupling)
    for i, oc in enumerate(focus_cols):
        orig_col_score.setdefault(int(oc), 0.0)
    return orig_col_score, focus_cols


# ----------------------------------------------------------------------------- public
def compute_coevolution(fasta_path, method="mi_apc", max_gap_frac=0.5,
                        weight_threshold=0.8, pseudocount=0.5, theta=0.2):
    """
    Returns dict:
      method, n_sequences, alignment_length, filtered_columns, neff,
      kept_cols (orig column indices), col_score (per filtered column),
      orig_col_score (dict orig_col -> score)
    """
    seqs = parse_msa(fasta_path)
    msa = encode_msa(seqs)
    msa_f, kept = filter_columns(msa, max_gap_frac=max_gap_frac)

    if method == "mi_apc":
        w = sequence_weights(msa_f, threshold=weight_threshold)
        neff = float(w.sum() / w.max())  # ~ effective sequence count proxy
        score = _mi_apc(msa_f, w, pseudocount=pseudocount)
        orig_col_score = {int(kept[i]): float(score[i]) for i in range(len(kept))}
        n_cols = int(msa_f.shape[1])
    elif method == "evcouplings":
        orig_col_score, focus_cols = _evcouplings_di(seqs, theta=theta)
        neff = None
        n_cols = len(focus_cols)
    else:
        raise ValueError(f"unknown method {method!r}; use 'mi_apc' or 'evcouplings'")

    return {
        "method": method,
        "n_sequences": len(seqs),
        "alignment_length": int(msa.shape[1]),
        "filtered_columns": n_cols,
        "neff": neff,
        "orig_col_score": {int(k): float(v) for k, v in orig_col_score.items()},
    }


def top_coupled_positions(coevo, col_to_pdb, fixed_already, n_top=42,
                          exclude_fixed=True):
    """
    Rank designable, positively-coupled positions and return the top-N PDB resnums.

    coevo         : output of compute_coevolution
    col_to_pdb    : ORIGINAL column -> PDB resnum (mapping.build_column_to_pdb)
    fixed_already : set of PDB resnums already fixed by conservation/contacts;
                    DCA's value is the positions these MISS, so exclude by default.
    """
    fixed_already = set(fixed_already)
    scored = []
    for orig_col, score in coevo["orig_col_score"].items():
        pdb = col_to_pdb.get(orig_col)
        if pdb is None or score <= 0:
            continue
        if exclude_fixed and pdb in fixed_already:
            continue
        scored.append((pdb, score))
    scored.sort(key=lambda x: x[1], reverse=True)
    top = scored[:n_top]
    return {
        "positions": sorted(p for p, _ in top),
        "score_cutoff": (top[-1][1] if top else 0.0),
        "ranked": [(p, round(s, 4)) for p, s in top],
    }
