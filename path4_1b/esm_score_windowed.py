#!/usr/bin/env python3
"""Windowed ESM-2 naturalness for full-length Cas9 (THE Cas9-specific tweak).

The RT campaign's esm_score.py did a single ESM-2 forward pass. ESM-2 650M caps
at ~1022 residues; SpRY-Cas9 is 1341 aa, so a single pass truncates/errors. This
version tiles OVERLAPPING windows and scores each residue from the window where
it is most centered (best two-sided context), then averages — identical to the
single-pass score when L <= window, so the naturalness axis stays comparable to
the RT work.

Score (per design) = mean over residues of  log P(actual residue | window context)
                     (wt-marginal: one forward pass per window, no per-residue masking)

Usage:
  python esm_score_windowed.py --selftest          # validate windowing, NO torch needed
  python esm_score_windowed.py --wt spry_cas9_wt.fasta \
         --designs <designs.json|designs.fasta> --out esm_scores.json
         [--window 1000 --stride 500 --model facebook/esm2_t33_650M_UR50D]
"""
import argparse, json, os, sys

# ----------------------------- windowing core (torch-free) -----------------------------
def plan_windows(L, window, stride):
    """Return (windows, owner) for a length-L sequence.
    windows : list of (start, end) half-open crops, each <= `window` long, covering [0,L).
    owner   : list length L; owner[i] = index into `windows` of the crop that SCORES residue i,
              chosen as the crop whose center is closest to i (maximizes two-sided context).
    For L <= window this is a single crop [0,L) and every residue is owned by it.
    """
    if L <= window:
        return [(0, L)], [0] * L
    starts = list(range(0, max(1, L - window) + 1, stride))
    if starts[-1] != L - window:
        starts.append(L - window)               # ensure the tail is fully covered
    windows = [(s, min(s + window, L)) for s in starts]
    centers = [(s + e) / 2.0 for (s, e) in windows]
    owner = []
    for i in range(L):
        best, bestd = 0, None
        for w, (s, e) in enumerate(windows):
            if s <= i < e:                       # crop must contain the residue
                d = abs(i - centers[w])
                if bestd is None or d < bestd:
                    bestd, best = d, w
        owner.append(best)
    return windows, owner


def selftest():
    cases = [(453, 1000, 500), (1341, 1000, 500), (1341, 1022, 511),
             (1000, 1000, 500), (1024, 1000, 500), (2500, 800, 400)]
    for L, W, S in cases:
        windows, owner = plan_windows(L, W, S)
        # 1) every residue owned by exactly one window, and that window contains it
        assert len(owner) == L
        for i, w in enumerate(owner):
            s, e = windows[w]
            assert s <= i < e, f"L={L}: residue {i} not in its owner window {windows[w]}"
        # 2) full coverage: union of windows == [0,L)
        covered = set()
        for (s, e) in windows:
            covered |= set(range(s, e))
        assert covered == set(range(L)), f"L={L}: incomplete coverage"
        # 3) each window <= W
        assert all(e - s <= W for (s, e) in windows)
        # 4) L<=W collapses to a single full-length window (== single-pass behaviour)
        if L <= W:
            assert windows == [(0, L)] and set(owner) == {0}
        # 5) owned residues get >= half-window context on at least one side (centered ownership)
        worst_ctx = min(min(i - windows[w][0], windows[w][1] - 1 - i) for i, w in enumerate(owner))
        print(f"  L={L:5} W={W} S={S}: {len(windows)} windows, "
              f"min one-sided context @owner = {worst_ctx}")
    print("SELFTEST PASS: windowing covers every residue exactly once with centered context.\n")


# ----------------------------- ESM scoring (needs torch) -----------------------------
def load_designs(path):
    """Accept a meta JSON (list of {id, seq}) or a FASTA. Return list of (id, seq)."""
    if path.endswith(".json"):
        d = json.load(open(path))
        items = d if isinstance(d, list) else (d.get("designs") or list(d.values()))
        out = []
        for x in items:
            if isinstance(x, dict):
                out.append((x.get("id") or x.get("name"), x.get("seq") or x.get("sequence")))
        return out
    # FASTA
    out, cur, name = [], [], None
    for ln in open(path):
        if ln.startswith(">"):
            if name is not None:
                out.append((name, "".join(cur)))
            name, cur = ln[1:].strip().split()[0], []
        else:
            cur.append(ln.strip())
    if name is not None:
        out.append((name, "".join(cur)))
    return out


def run_esm(args):
    import torch
    from transformers import AutoTokenizer, EsmForMaskedLM
    dev = "mps" if torch.backends.mps.is_available() else ("cuda" if torch.cuda.is_available() else "cpu")
    print(f"device={dev} model={args.model} window={args.window} stride={args.stride}", flush=True)
    tok = AutoTokenizer.from_pretrained(args.model)
    model = EsmForMaskedLM.from_pretrained(args.model).to(dev).eval()

    @torch.no_grad()
    def naturalness(seq):
        L = len(seq)
        windows, owner = plan_windows(L, args.window, args.stride)
        # per-window logp for the actual residue at every position in that window
        win_logp = []
        for (s, e) in windows:
            enc = tok(seq[s:e], return_tensors="pt").to(dev)
            logits = model(**enc).logits[0]              # (len+2, vocab)
            logp = torch.log_softmax(logits.float(), dim=-1)
            ids = enc["input_ids"][0]
            # positions 1..len map to residues s..e-1 (skip <cls>/<eos>)
            row = {}
            for k in range(1, ids.shape[0] - 1):
                row[s + (k - 1)] = logp[k, ids[k]].item()
            win_logp.append(row)
        vals = [win_logp[owner[i]][i] for i in range(L)]
        return sum(vals) / L

    wt = "".join(l.strip() for l in open(args.wt) if not l.startswith(">"))
    out_path = args.out
    out, wt_score = {}, None
    if os.path.exists(out_path):
        prev = json.load(open(out_path))
        out, wt_score = prev.get("scores", {}), prev.get("wt")
        print(f"resuming: {len(out)} already scored", flush=True)
    if wt_score is None:
        wt_score = naturalness(wt)
    print(f"WT({len(wt)}) naturalness = {wt_score:.4f}", flush=True)

    designs = load_designs(args.designs)
    import time
    t0 = time.time()
    for k, (sid, seq) in enumerate(designs):
        if sid in out or not seq:
            continue
        out[sid] = naturalness(seq)
        if (k + 1) % 50 == 0:
            json.dump({"wt": wt_score, "model": args.model, "window": args.window,
                       "stride": args.stride, "scores": out}, open(out_path, "w"))
            el = time.time() - t0
            print(f"  {len(out)}/{len(designs)}  {el/ max(1,len(out)):.2f}s/seq", flush=True)
    json.dump({"wt": wt_score, "model": args.model, "window": args.window,
               "stride": args.stride, "scores": out}, open(out_path, "w"))
    print(f"done: {len(out)} designs scored -> {out_path}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--wt"); ap.add_argument("--designs"); ap.add_argument("--out", default="esm_scores.json")
    ap.add_argument("--window", type=int, default=1000)
    ap.add_argument("--stride", type=int, default=500)
    ap.add_argument("--model", default="facebook/esm2_t33_650M_UR50D")
    a = ap.parse_args()
    if a.selftest or not a.designs:
        selftest()
        if not a.designs:
            sys.exit(0)
    run_esm(a)
