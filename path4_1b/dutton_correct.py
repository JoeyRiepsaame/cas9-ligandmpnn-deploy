#!/usr/bin/env python3
"""Dutton composition-bias correction of MPNN log-odds (post-generation step 2).

The principled replacement for the dropped E-penalty (RT-4.1b). MPNN has a
per-amino-acid compositional bias; Dutton de-biases the score by subtracting that
background so ranking reflects position-specific fit, not the model's global AA
preference.

  ref[a]          = mean over ALL designs & designable positions of log_probs[:,i,a]
                    (the model's background log-prob for amino acid a)
  raw(design)     = mean_i log P(actual residue_i)                 over designable i
  corrected(design)= mean_i ( log P(residue_i) - ref[residue_i] )  over designable i

Reads the LigandMPNN stats .pt tree, joins to unique_meta.json by sequence, writes
dutton_scores.json {ref, designs:[{id, raw_mpnn, mpnn_corr}]}.

Usage: python dutton_correct.py --root <outputs> --meta unique_meta.json --out dutton_scores.json
"""
import argparse, glob, json, os
import numpy as np, torch

ALPHA = "ACDEFGHIKLMNPQRSTVWYX"; A = {a: i for i, a in enumerate(ALPHA)}


def main(args):
    pts = sorted(glob.glob(os.path.join(args.root, "**", "stats", "*.pt"), recursive=True))
    # pass 1: accumulate per-AA background over designable positions
    bg_sum = np.zeros(21); bg_n = 0
    designable = None
    for pt in pts:
        d = torch.load(pt, map_location="cpu", weights_only=False)
        lp = np.asarray(d["log_probs"])                  # (B,L,21)
        cm = np.asarray(d["chain_mask"]).astype(bool)
        designable = cm
        idx = np.where(cm)[0]
        bg_sum += lp[:, idx, :].reshape(-1, 21).sum(0)
        bg_n += lp.shape[0] * len(idx)
    ref = bg_sum / max(1, bg_n)                          # (21,) background log-prob per AA
    print("Dutton background ref (per-AA, designable positions):")
    print("  " + "  ".join(f"{ALPHA[i]}:{ref[i]:.2f}" for i in range(20)))

    # pass 2: per-design raw + corrected, keyed by sequence
    by_seq = {}
    for pt in pts:
        d = torch.load(pt, map_location="cpu", weights_only=False)
        gs = np.asarray(d["generated_sequences"]); lp = np.asarray(d["log_probs"])
        cm = np.asarray(d["chain_mask"]).astype(bool); idx = np.where(cm)[0]
        for b in range(gs.shape[0]):
            seq = "".join(ALPHA[c] for c in gs[b])
            raw = float(np.mean([lp[b, i, gs[b, i]] for i in idx]))
            corr = float(np.mean([lp[b, i, gs[b, i]] - ref[gs[b, i]] for i in idx]))
            # keep best (highest corrected) per unique sequence
            if seq not in by_seq or corr > by_seq[seq][1]:
                by_seq[seq] = (round(raw, 4), round(corr, 4))

    meta = json.load(open(args.meta))
    out = []
    miss = 0
    for r in meta:
        v = by_seq.get(r["seq"])
        if v is None:
            miss += 1; continue
        out.append(dict(id=r["id"], tier=r["tier"], raw_mpnn=v[0], mpnn_corr=v[1]))
    print(f"scored {len(out)}/{len(meta)} unique designs (join misses: {miss})")
    json.dump({"ref": {ALPHA[i]: round(float(ref[i]), 4) for i in range(20)},
               "designs": out}, open(args.out, "w"))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--meta", default="unique_meta.json")
    ap.add_argument("--out", default="dutton_scores.json")
    main(ap.parse_args())
