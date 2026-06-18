"""COMPSS Pareto ranking + synthesis shortlist (post-generation step 3).

Combines the orthogonal quality axes the RT campaign validated:
  ESM naturalness (windowed)  ⟂  Dutton-corrected MPNN  ( +optional design-LD score)
and %identity-to-WT for the patent-novelty gate. Designs on the ESM×Dutton Pareto
frontier that also pass the identity ceiling form the synthesis shortlist.

Inputs (by design id):
  unique_meta.json      (id, seq, tier)
  esm_scores.json       (from esm_score_windowed.py;  {wt, scores:{id:val}})
  dutton_scores.json    (from dutton_correct.py)
  [ld_scoring_table.json optional: design_ld_scores]
WT reference: spry_cas9_wt.fasta. Writes pareto_ranked.csv + synthesis_shortlist.json.
Sequences go to the shortlist FILE only (audit before synthesis), never printed.

Usage: python pareto_rank.py [--id-ceiling 70 --top 24]
"""
import argparse, csv, json, os, statistics as st
from collections import defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))


def pct_identity(wt, design):
    from Bio import Align
    al = Align.PairwiseAligner(); al.mode = "global"
    al.open_gap_score = -10; al.extend_gap_score = -0.5
    al.match_score = 2; al.mismatch_score = -1
    aln = al.align(wt, design)[0]
    ident = aligned = 0
    for (ws, we), (ds, de) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(we - ws):
            aligned += 1; ident += (wt[ws + k] == design[ds + k])
    return round(100.0 * ident / aligned, 1)


def model_of(best_src):
    """model type from extract_seqs best_src ('<tier>/<run>#<b>'): sol/prot/lig -> name."""
    p = best_src.split("/")[1].split("#")[0].split("_")[0]
    return {"sol": "Soluble", "prot": "Protein", "lig": "Ligand"}.get(p, p)


def pareto_front(rows, ax1, ax2):
    front = []
    for r in rows:
        if not any((o[ax1] >= r[ax1] and o[ax2] >= r[ax2] and
                    (o[ax1] > r[ax1] or o[ax2] > r[ax2])) for o in rows):
            front.append(r)
    return front


def main(a):
    meta = {r["id"]: r for r in json.load(open(os.path.join(HERE, a.meta)))}
    esm = json.load(open(os.path.join(HERE, a.esm)))
    esm_s, wt_esm = esm["scores"], esm["wt"]
    dutton = {d["id"]: d for d in json.load(open(os.path.join(HERE, a.dutton)))["designs"]}
    ld = {}
    if os.path.exists(os.path.join(HERE, a.ld)):
        ld = json.load(open(os.path.join(HERE, a.ld))).get("design_ld_scores", {})
    wt = "".join(l.strip() for l in open(os.path.join(HERE, "spry_cas9_wt.fasta")) if not l.startswith(">"))

    rows = []
    for sid, m in meta.items():
        if sid not in esm_s or sid not in dutton:
            continue
        rows.append(dict(id=sid, tier=m["tier"], model=model_of(m["best_src"]), esm=esm_s[sid],
                         mpnn_corr=dutton[sid]["mpnn_corr"], mpnn_raw=dutton[sid]["raw_mpnn"],
                         ld=ld.get(sid), pct_identity=pct_identity(wt, m["seq"]), seq=m["seq"]))
    print(f"merged {len(rows)} designs (WT esm={wt_esm:.4f})")
    if not rows:
        print("no designs joined — check esm/dutton inputs"); return

    # FAIRNESS FIX (per-model z): the Dutton/MPNN score is NOT comparable across model
    # types — LigandMPNN conditions on atom context, inflating its self-confidence vs
    # Soluble/Protein. z-normalise mpnn_corr WITHIN each model type so all three compete
    # fairly. ESM is model-agnostic, so its raw value is kept for the Pareto axis.
    by_model = defaultdict(list)
    for r in rows:
        by_model[r["model"]].append(r["mpnn_corr"])
    mstat = {m: (st.mean(v), st.pstdev(v) or 1) for m, v in by_model.items()}
    for r in rows:
        mu, sd = mstat[r["model"]]
        r["z_mpnn"] = (r["mpnn_corr"] - mu) / sd
    print("per-model mpnn_corr mean (pre-norm bias): "
          + ", ".join(f"{m}={mstat[m][0]:.3f}(n={len(by_model[m])})" for m in sorted(mstat)))

    # Pareto + combined z use the per-model-normalised MPNN axis vs model-agnostic ESM
    front = pareto_front(rows, "esm", "z_mpnn")
    me, se = st.mean(r["esm"] for r in rows), st.pstdev(r["esm"] for r in rows) or 1
    for r in rows:
        r["z"] = (r["esm"] - me) / se + r["z_mpnn"]
    rows.sort(key=lambda r: -r["z"])

    ids = sorted(r["pct_identity"] for r in rows)
    print(f"Pareto frontier (esm x Dutton): {len(front)} | "
          f"%id min/median/max = {ids[0]:.1f}/{ids[len(ids)//2]:.1f}/{ids[-1]:.1f} | "
          f"<{a.id_ceiling}%: {sum(1 for x in ids if x < a.id_ceiling)}")
    for grp, key in (("tier", "tier"), ("model", "model")):
        pt = defaultdict(list)
        for r in rows:
            pt[r[key]].append(r)
        print(f"per-{grp}:")
        for t in sorted(pt):
            v = pt[t]
            print(f"  {t:14} n={len(v):4} esm={st.mean(x['esm'] for x in v):.3f} "
                  f"mpnn_corr={st.mean(x['mpnn_corr'] for x in v):.3f} "
                  f"%id={st.mean(x['pct_identity'] for x in v):.1f}")

    with open(os.path.join(HERE, "pareto_ranked.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "tier", "model", "esm", "mpnn_corr", "z_mpnn", "mpnn_raw",
                    "ld", "pct_identity", "z"])
        for r in rows:
            w.writerow([r["id"], r["tier"], r["model"], r["esm"], r["mpnn_corr"],
                        round(r["z_mpnn"], 3), r["mpnn_raw"], r["ld"], r["pct_identity"],
                        round(r["z"], 3)])
    front_ids = {r["id"] for r in front}
    shortlist = [r for r in rows if r["id"] in front_ids and r["pct_identity"] < a.id_ceiling][:a.top]
    from collections import Counter
    print(f"\nPareto frontier {len(front)} | shortlist {len(shortlist)} (Pareto ∩ <{a.id_ceiling}% id)")
    print("  shortlist by MODEL:", dict(Counter(r["model"] for r in shortlist)),
          "| by tier:", dict(Counter(r["tier"] for r in shortlist)))
    json.dump({"wt_esm": wt_esm, "id_ceiling": a.id_ceiling, "pareto_frontier": len(front),
               "ranking": "ESM (model-agnostic) x per-model-z Dutton-MPNN",
               "shortlist": [{k: r[k] for k in ("id", "tier", "model", "esm", "mpnn_corr",
                              "z_mpnn", "pct_identity", "seq")} for r in shortlist]},
              open(os.path.join(HERE, "synthesis_shortlist.json"), "w"), indent=2)
    print(f"wrote pareto_ranked.csv + synthesis_shortlist.json")
    print("NEXT: audit shortlist sequences, then AF3 triage (with calibrators) + wet-lab PAM panel.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", default="unique_meta.json")
    ap.add_argument("--esm", default="esm_scores.json")
    ap.add_argument("--dutton", default="dutton_scores.json")
    ap.add_argument("--ld", default="ld_scoring_table.json")
    ap.add_argument("--id-ceiling", type=float, default=70.0, dest="id_ceiling")
    ap.add_argument("--top", type=int, default=24)
    main(ap.parse_args())
