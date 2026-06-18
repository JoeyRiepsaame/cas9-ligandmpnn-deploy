#!/usr/bin/env python3
"""Extract + dedup the generated SpRY-Cas9 designs (post-generation step 1).

Reads every LigandMPNN stats .pt under the output tree (each holds
generated_sequences (B,1341), log_probs (B,1341,21), native_sequence, chain_mask)
and the matching seqs/*.fa header (overall_confidence per design). Decodes each
design to a full chain-A sequence, dedups, and writes:

  unique_meta.json   [{id, seq, tier, n_copies, best_overall_confidence,
                       best_src, raw_mpnn, aa_after_omitC_check}]
  unique.fasta

raw_mpnn = mean log P(actual residue) over DESIGNABLE positions (chain_mask) —
the uncorrected MPNN score; dutton_correct.py de-biases it. Sequences are written
to file (audited downstream), never printed.

Usage: python extract_seqs.py --root <outputs_dir> [--out-meta unique_meta.json --out-fa unique.fasta]
"""
import argparse, glob, json, os, re
import numpy as np, torch

ALPHA = "ACDEFGHIKLMNPQRSTVWYX"

def parse_fa_conf(fa_path):
    """map design index (1-based id=) -> overall_confidence from the .fa headers."""
    conf = {}
    if not os.path.exists(fa_path):
        return conf
    for ln in open(fa_path):
        if ln.startswith(">") and "id=" in ln:
            mid = re.search(r"id=(\d+)", ln); mc = re.search(r"overall_confidence=([0-9.]+)", ln)
            if mid and mc:
                conf[int(mid.group(1))] = float(mc.group(1))
    return conf


def main(a):
    pts = sorted(glob.glob(os.path.join(a.root, "**", "stats", "*.pt"), recursive=True))
    print(f"found {len(pts)} stats .pt files under {a.root}")
    uniq = {}            # seq -> meta
    n_designs = 0
    for pt in pts:
        rel = os.path.relpath(pt, a.root)
        tier = rel.split(os.sep)[0]
        src = os.path.dirname(os.path.dirname(rel))      # <tier>/<run>
        d = torch.load(pt, map_location="cpu", weights_only=False)
        gs = np.asarray(d["generated_sequences"])         # (B,1341)
        lp = np.asarray(d["log_probs"])                   # (B,1341,21)
        cm = np.asarray(d["chain_mask"]).astype(bool)     # (1341,) designable
        fa = pt.replace(os.sep + "stats" + os.sep, os.sep + "seqs" + os.sep).rsplit(".pt", 1)[0] + ".fa"
        conf = parse_fa_conf(fa)
        for b in range(gs.shape[0]):
            n_designs += 1
            seq = "".join(ALPHA[c] for c in gs[b])
            # raw MPNN score: mean log P(actual residue) over designable positions
            idx = np.where(cm)[0]
            raw = float(np.mean([lp[b, i, gs[b, i]] for i in idx]))
            oc = conf.get(b + 1)
            m = uniq.get(seq)
            if m is None:
                uniq[seq] = dict(seq=seq, tier=tier, n_copies=1,
                                 best_overall_confidence=oc, best_src=f"{src}#{b}",
                                 raw_mpnn=round(raw, 4))
            else:
                m["n_copies"] += 1
                if oc is not None and (m["best_overall_confidence"] is None
                                       or oc > m["best_overall_confidence"]):
                    m["best_overall_confidence"] = oc; m["best_src"] = f"{src}#{b}"
                if raw > m["raw_mpnn"]:
                    m["raw_mpnn"] = round(raw, 4)
    rows = list(uniq.values())
    for i, r in enumerate(rows):
        r["id"] = f"d{i:05d}"
    # quick omit-C integrity: every design should have C only where WT-fixed (574)
    bad_c = sum(1 for r in rows if r["seq"].count("C") > 1)
    print(f"designs: {n_designs} total, {len(rows)} unique "
          f"({100*len(rows)/max(1,n_designs):.0f}% unique)")
    print(f"per-tier unique: " + ", ".join(
        f"{t}={sum(1 for r in rows if r['tier']==t)}" for t in sorted(set(r['tier'] for r in rows))))
    print(f"designs with >1 Cys (should be 0; only fixed C574 allowed): {bad_c}")
    json.dump(rows, open(a.out_meta, "w"))
    with open(a.out_fa, "w") as fh:
        for r in rows:
            fh.write(f">{r['id']} tier={r['tier']} oc={r['best_overall_confidence']}\n")
            for i in range(0, len(r["seq"]), 60):
                fh.write(r["seq"][i:i+60] + "\n")
    print(f"wrote {a.out_meta} + {a.out_fa}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--out-meta", default="unique_meta.json", dest="out_meta")
    ap.add_argument("--out-fa", default="unique.fasta", dest="out_fa")
    main(ap.parse_args())
