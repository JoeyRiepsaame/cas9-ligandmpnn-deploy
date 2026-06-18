#!/usr/bin/env python3
"""Design-ensemble LD / coevolution for SpRY-Cas9 (gap C2/C3).

DISTINCT from the MSA-DCA in constraints/: that measures EVOLUTIONARY coupling
across natural Cas9 orthologs. This measures coupling WITHIN the MPNN design
ensemble — which free positions co-vary across the generated designs. High
design-LD pairs are either real epistasis the model learned, or constraint-
induced artifacts (the RT campaign's C90-184 was the latter). We also build a
joint-frequency table to SCORE each design for internal LD-compatibility.

Cas9-scale tweak: the blind-spot detector's compute_mi_matrix is pure-Python
(O(P²·N) loops) — fine for RT's 60 designs, far too slow for 10K designs × ~600
variable positions. This is a vectorised MI+APC (one-hot block matmul).

Usage:
  python ld_analysis.py --selftest                       # validate MI+APC+scoring, no inputs
  python ld_analysis.py --designs <fasta|meta.json> \
        --map wt_resnum_to_index.json \
        --out-ld ld_matrix.json --out-score ld_scoring_table.json \
        [--min-maf 0.02 --top 100 --max-pos 1200]
"""
import argparse, json, os, sys
import numpy as np

ALPHA = "ACDEFGHIKLMNPQRSTVWY"
A2I = {a: i for i, a in enumerate(ALPHA)}
NA = 20


def load_designs(path):
    if path.endswith(".json"):
        d = json.load(open(path))
        items = d if isinstance(d, list) else (d.get("designs") or list(d.values()))
        return [(x.get("id") or x.get("name"), x.get("seq") or x.get("sequence"))
                for x in items if isinstance(x, dict)]
    out, cur, name = [], [], None
    for ln in open(path):
        if ln.startswith(">"):
            if name is not None:
                out.append((name, "".join(cur)))
            name, cur = ln[1:].split()[0], []
        else:
            cur.append(ln.strip())
    if name is not None:
        out.append((name, "".join(cur)))
    return out


def encode(seqs):
    """(N,L) int8; AA 0-19, anything else (incl. gap) = -1."""
    N, L = len(seqs), len(seqs[0])
    X = np.full((N, L), -1, dtype=np.int8)
    for n, s in enumerate(seqs):
        for j, c in enumerate(s):
            X[n, j] = A2I.get(c, -1)
    return X


def variable_positions(X, min_maf):
    """columns whose 2nd-most-common AA freq >= min_maf (i.e. genuinely variable)."""
    N = X.shape[0]
    out = []
    for c in range(X.shape[1]):
        col = X[:, c]
        col = col[col >= 0]
        if col.size == 0:
            continue
        counts = np.bincount(col, minlength=NA)
        srt = np.sort(counts)[::-1]
        if srt.size >= 2 and srt[1] / N >= min_maf:
            out.append(c)
    return out


def mi_apc_vectorized(X, positions):
    """Vectorised MI + APC over `positions`. Returns (mi_apc, f_single dict, onehot info)."""
    N = X.shape[0]
    P = len(positions)
    # one-hot block matrix Of: (N, P*NA)
    Of = np.zeros((N, P * NA), dtype=np.float32)
    for pi, c in enumerate(positions):
        col = X[:, c]
        valid = col >= 0
        cols = pi * NA + col[valid].astype(np.intp)       # intp: avoid int8 overflow at scale
        Of[np.where(valid)[0], cols] = 1.0
    # macOS Accelerate emits spurious FP warnings on float32 matmul; the result is
    # correct (asserted finite below), so silence only the cosmetic warning here.
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        G = (Of.T @ Of) / N                              # (P*NA, P*NA) all joint freqs
    assert np.isfinite(G).all(), "non-finite joint-frequency matrix (real numerical error)"
    assert G.min() >= -1e-6 and G.max() <= 1.0 + 1e-6, "joint freqs out of [0,1] — bug"
    # single freqs = diagonal of each (i,i) block
    f1 = np.zeros((P, NA), dtype=np.float64)
    for i in range(P):
        f1[i] = np.diag(G[i*NA:(i+1)*NA, i*NA:(i+1)*NA])
    mi = np.zeros((P, P), dtype=np.float64)
    for i in range(P):
        fi = f1[i]
        for j in range(i + 1, P):
            f2 = G[i*NA:(i+1)*NA, j*NA:(j+1)*NA].astype(np.float64)   # (NA,NA) joint
            outer = np.outer(fi, f1[j])
            mask = (f2 > 0) & (outer > 0)
            val = float((f2[mask] * np.log2(f2[mask] / outer[mask])).sum())
            mi[i, j] = mi[j, i] = val
    # APC
    row_means = mi.mean(axis=1)
    gmean = mi.mean()
    mi_apc = mi.copy()
    if gmean > 0:
        apc = np.outer(row_means, row_means) / gmean
        mi_apc = mi - apc
        np.fill_diagonal(mi_apc, 0.0)
    return mi_apc, f1, G


def top_pairs(mi_apc, positions, idx2label, k, anti_k):
    P = len(positions)
    pairs = [(i, j, mi_apc[i, j]) for i in range(P) for j in range(i + 1, P)]
    pairs.sort(key=lambda t: -t[2])
    def fmt(i, j, v):
        return {"pos_i": idx2label(positions[i]), "pos_j": idx2label(positions[j]),
                "mi_apc": round(float(v), 4)}
    return ([fmt(*p) for p in pairs[:k]], [fmt(*p) for p in pairs[-anti_k:][::-1]])


def build_scoring_table(X, positions, top, idx2label, label2idx):
    """For each top pair, store the joint AA-frequency table P(a,b) for design scoring."""
    N = X.shape[0]
    table = []
    for pr in top:
        i = label2idx[pr["pos_i"]]; j = label2idx[pr["pos_j"]]
        ci, cj = X[:, i], X[:, j]
        joint = {}
        for a in range(NA):
            for b in range(NA):
                cnt = int(np.sum((ci == a) & (cj == b)))
                if cnt:
                    joint[f"{ALPHA[a]}{ALPHA[b]}"] = cnt / N
        table.append({"pos_i": pr["pos_i"], "pos_j": pr["pos_j"],
                      "i_idx": i, "j_idx": j, "mi_apc": pr["mi_apc"], "joint": joint})
    return table


def score_designs(X, table, ids, eps=1e-4):
    """LD-compatibility per design = mean log2(P(aa_i,aa_j)+eps) over scoring pairs."""
    scores = {}
    for n in range(X.shape[0]):
        tot, m = 0.0, 0
        for t in table:
            a = X[n, t["i_idx"]]; b = X[n, t["j_idx"]]
            if a < 0 or b < 0:
                continue
            p = t["joint"].get(f"{ALPHA[a]}{ALPHA[b]}", 0.0)
            tot += np.log2(p + eps); m += 1
        scores[ids[n]] = round(tot / m, 4) if m else None
    return scores


# ------------------------------------------------------------------ selftest
def selftest():
    rng_seq = []
    # synthetic ensemble: 800 designs, 6 positions. Plant a STRONG coupling between
    # pos1 and pos2 (always 'KE' or 'RD'), independence elsewhere; pos0 invariant.
    N = 800
    for n in range(N):
        s = list("A______")
        # coupled block: 50/50 between the two allowed combos
        if n % 2 == 0:
            s[1], s[2] = "K", "E"
        else:
            s[1], s[2] = "R", "D"
        # independent variable positions 3,4 (no coupling)
        s[3] = "L" if (n // 2) % 2 == 0 else "V"
        s[4] = "S" if (n // 3) % 2 == 0 else "T"
        s[5] = "G"                                  # invariant
        rng_seq.append("".join(c if c != "_" else "A" for c in s))
    X = encode(rng_seq)
    var = variable_positions(X, min_maf=0.05)
    assert set(var) == {1, 2, 3, 4}, f"variable positions wrong: {var}"
    mi_apc, f1, _ = mi_apc_vectorized(X, var)
    lbl = lambda p: f"p{p}"
    top, anti = top_pairs(mi_apc, var, lbl, k=1, anti_k=1)
    # the planted 1-2 coupling must be the #1 pair
    assert {top[0]["pos_i"], top[0]["pos_j"]} == {"p1", "p2"}, f"top pair not 1-2: {top[0]}"
    l2i = {lbl(p): p for p in var}
    table = build_scoring_table(X, var, top, lbl, l2i)
    # a COMPATIBLE design (KE) should score higher than an INCOMPATIBLE one (KD, never seen)
    good = encode(["AKE LSG".replace(" ", "")])
    bad = encode(["AKD LSG".replace(" ", "")])
    sg = score_designs(good, table, ["good"])["good"]
    sb = score_designs(bad, table, ["bad"])["bad"]
    assert sg > sb, f"LD scoring failed to penalise incompatible design: good={sg} bad={sb}"
    # cross-check vectorised MI vs the pure-python detector on this small case
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from mpnn_blind_spot_detector import compute_mi_matrix
    _, mi_ref, _ = compute_mi_matrix(rng_seq, var)
    # compare raw MI ordering of the 1-2 pair being max (APC differs in detail; check MI signal)
    i, j = var.index(1), var.index(2)
    assert mi_ref[i, j] == max(mi_ref[i, k] for k in range(len(var)) if k != i) , "ref MI disagrees"
    print(f"SELFTEST PASS: planted 1-2 coupling recovered as #1 (mi_apc={top[0]['mi_apc']}); "
          f"LD score good={sg} > bad={sb}; vectorised MI matches detector on the coupled pair.\n")


def run(args):
    designs = [(i, s) for i, s in load_designs(args.designs) if s]
    ids = [d[0] for d in designs]
    seqs = [d[1] for d in designs]
    L = len(seqs[0])
    assert all(len(s) == L for s in seqs), "designs not equal length"
    X = encode(seqs)
    print(f"ensemble: {len(seqs)} designs x {L} positions")
    var = variable_positions(X, args.min_maf)
    print(f"variable positions (MAF>={args.min_maf}): {len(var)}")
    if len(var) > args.max_pos:
        print(f"  capping to top-{args.max_pos} most-variable (memory). "
              f"P*20={args.max_pos*NA} -> G ~ {(args.max_pos*NA)**2*4/1e9:.1f}GB float32")
        # keep highest minor-allele-freq positions
        maf = []
        for c in var:
            col = X[:, c]; col = col[col >= 0]
            srt = np.sort(np.bincount(col, minlength=NA))[::-1]
            maf.append((srt[1] / len(seqs), c))
        var = [c for _, c in sorted(maf, reverse=True)[:args.max_pos]]
    # map position index -> author resnum if a map is given
    idx2author = None
    if args.map and os.path.exists(args.map):
        r2i = {int(k): v for k, v in json.load(open(args.map)).items()}
        idx2author = {v: k for k, v in r2i.items()}
    lbl = (lambda p: f"A{idx2author[p]}" if idx2author and p in idx2author else f"i{p}")
    mi_apc, f1, _ = mi_apc_vectorized(X, var)
    top, anti = top_pairs(mi_apc, var, lbl, args.top, args.top // 2)
    l2i = {lbl(p): p for p in var}
    table = build_scoring_table(X, var, top, lbl, l2i)
    scores = score_designs(X, table, ids)
    iu = np.triu_indices(len(var), 1)
    json.dump({"n_designs": len(seqs), "n_variable_positions": len(var),
               "variable_positions": [lbl(p) for p in var],
               "mi_apc_mean": round(float(mi_apc[iu].mean()), 5),
               "top_pairs": top, "anti_pairs": anti},
              open(args.out_ld, "w"), indent=2)
    json.dump({"description": "LD scoring table; design score = mean log2(P(aa_i,aa_j)+eps) over pairs",
               "n_designs": len(seqs), "pairs": table, "design_ld_scores": scores},
              open(args.out_score, "w"), indent=2)
    print(f"wrote {args.out_ld} ({len(top)} top pairs) + {args.out_score} (LD scores for {len(scores)} designs)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--designs"); ap.add_argument("--map")
    ap.add_argument("--out-ld", default="ld_matrix.json", dest="out_ld")
    ap.add_argument("--out-score", default="ld_scoring_table.json", dest="out_score")
    ap.add_argument("--min-maf", type=float, default=0.02, dest="min_maf")
    ap.add_argument("--top", type=int, default=100)
    ap.add_argument("--max-pos", type=int, default=1200, dest="max_pos")
    a = ap.parse_args()
    if a.selftest or not a.designs:
        selftest()
        if not a.designs:
            sys.exit(0)
    run(a)
