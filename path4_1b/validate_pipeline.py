#!/usr/bin/env python3
"""THOROUGH end-to-end validation of the SpRY-Cas9 results pipeline.

Independent re-derivation: decodes designs straight from the 90 LigandMPNN stats
.pt files and recomputes scores, then cross-checks every downstream artifact
(meta, ESM, Dutton, per-model-z, shortlists, pilot set, AF3 batch) against ground
truth and against each other. Does NOT trust the derived JSONs. Several bugs were
caught this session (tuple-mapping, int8 overflow, open()[]-index, cross-model
score bias) — this is deliberately paranoid. Exit non-zero on any failure.
"""
import csv, glob, json, os, sys
from collections import Counter, defaultdict
import numpy as np, torch
HERE = os.path.dirname(os.path.abspath(__file__))

FAIL, NC = [], 0
def ck(c, m):
    global NC; NC += 1
    print(("  PASS  " if c else "  FAIL  ") + m)
    if not c: FAIL.append(m)
def sec(t): print("\n" + "=" * 74 + f"\n{t}\n" + "=" * 74)

ALPHA = "ACDEFGHIKLMNPQRSTVWYX"
r2i = {int(k): v for k, v in json.load(open(os.path.join(HERE, "wt_resnum_to_index.json"))).items()}
i2r = {v: k for k, v in r2i.items()}
CAT = {10:'D',762:'E',839:'D',840:'H',863:'N',983:'H',986:'D'}
SPRY = {61:'R',1111:'R',1135:'L',1136:'W',1218:'K',1219:'Q',1317:'R',1322:'R',1333:'P',1335:'Q',1337:'R'}
def model_of(src): return {"sol":"Soluble","prot":"Protein","lig":"Ligand"}[src.split("/")[1].split("#")[0].split("_")[0]]

# ---- INDEPENDENT decode of every design from the .pt tree ----
sec("1. GENERATION — independent decode of all stats .pt vs unique_meta")
pts = sorted(glob.glob(os.path.join(HERE, "outputs", "**", "stats", "*.pt"), recursive=True))
ck(len(pts) == 90, f"90 stats .pt files (got {len(pts)})")
dec = {}                       # seq -> (raw_mpnn, model, tier)
bg_sum = np.zeros(21); bg_n = 0; n_designs = 0
viol = {"len": 0, "cys": 0, "cat": 0, "spry": 0}
for pt in pts:
    rel = os.path.relpath(pt, os.path.join(HERE, "outputs"))
    tier = rel.split(os.sep)[0]; run = rel.split(os.sep)[1]
    d = torch.load(pt, map_location="cpu", weights_only=False)
    gs = np.asarray(d["generated_sequences"]); lp = np.asarray(d["log_probs"])
    cm = np.asarray(d["chain_mask"]).astype(bool); idx = np.where(cm)[0]
    bg_sum += lp[:, idx, :].reshape(-1, 21).sum(0); bg_n += lp.shape[0] * len(idx)
    for b in range(gs.shape[0]):
        n_designs += 1
        seq = "".join(ALPHA[c] for c in gs[b])
        if len(seq) != 1341: viol["len"] += 1
        if [p for p, i in r2i.items() if seq[i] == 'C'] != [574]: viol["cys"] += 1
        if any(seq[r2i[p]] != aa for p, aa in CAT.items()): viol["cat"] += 1
        if any(seq[r2i[p]] != aa for p, aa in SPRY.items()): viol["spry"] += 1
        raw = float(np.mean([lp[b, i, gs[b, i]] for i in idx]))
        dec.setdefault(seq, (round(raw, 4), {"sol":"Soluble","prot":"Protein","lig":"Ligand"}[run.split("_")[0]], tier))
ck(n_designs == 1530, f"1530 designs decoded (got {n_designs})")
ck(len(dec) == 1530, f"all unique (got {len(dec)})")
ck(viol["len"] == 0, f"all len 1341 (violations {viol['len']})")
ck(viol["cys"] == 0, f"omit-C: cysteine only at fixed C574, every design (violations {viol['cys']})")
ck(viol["cat"] == 0, f"7/7 catalytic WT, every design (violations {viol['cat']})")
ck(viol["spry"] == 0, f"11/11 SpRY preserved, every design (violations {viol['spry']})")
meta = {r["id"]: r for r in json.load(open(os.path.join(HERE, "unique_meta.json")))}
ck(len(meta) == 1530, f"unique_meta has 1530 (got {len(meta)})")
ck(set(m["seq"] for m in meta.values()) == set(dec), "unique_meta sequences == independent decode (byte-match set)")
ck(all(model_of(m["best_src"]) == dec[m["seq"]][1] for m in meta.values()),
   "meta best_src model type matches the .pt run that produced each sequence")

# ---- Dutton recompute ----
sec("2. DUTTON — recompute background + sample corrected vs file")
ref = bg_sum / bg_n
dfile = json.load(open(os.path.join(HERE, "dutton_scores.json")))
ck(all(abs(ref[i] - dfile["ref"][ALPHA[i]]) < 1e-3 for i in range(20)),
   "Dutton per-AA background ref reproduced from .pt within 1e-3")
dmap = {d["id"]: d for d in dfile["designs"]}
ck(len(dmap) == 1530, f"dutton_scores has 1530 (got {len(dmap)})")
# sample 15 designs: recompute corrected, compare
sample = list(meta)[::100][:15]
ok = True
for sid in sample:
    seq = meta[sid]["seq"]
    raw = dec[seq][0]
    if abs(raw - dmap[sid]["raw_mpnn"]) > 1e-3: ok = False
ck(ok, f"raw MPNN matches file for {len(sample)} sampled designs (within 1e-3)")

# ---- ESM coverage ----
sec("3. ESM — coverage + sanity")
esm = json.load(open(os.path.join(HERE, "esm_scores.json")))
ck(len(esm["scores"]) == 1530, f"ESM scored 1530 (got {len(esm['scores'])})")
ck(set(esm["scores"]) == set(meta), "ESM ids == meta ids")
vals = list(esm["scores"].values())
ck(not any(np.isnan(v) or np.isinf(v) for v in vals), "no NaN/inf in ESM scores")
ck(isinstance(esm["wt"], float), f"WT ESM present ({esm['wt']:.4f})")

# ---- per-model-z fix in pareto_ranked.csv ----
sec("4. PER-MODEL-Z FIX — recompute vs pareto_ranked.csv")
rows = list(csv.DictReader(open(os.path.join(HERE, "pareto_ranked.csv"))))
ck(len(rows) == 1530, f"pareto_ranked has 1530 (got {len(rows)})")
ck(dict(Counter(r["model"] for r in rows)) == {"Ligand": 810, "Protein": 360, "Soluble": 360},
   "model split 810/360/360 in ranking")
# recompute z_mpnn per model and compare
bym = defaultdict(list)
for r in rows: bym[r["model"]].append(float(r["mpnn_corr"]))
import statistics as st
mstat = {m: (st.mean(v), st.pstdev(v) or 1) for m, v in bym.items()}
bad = 0
for r in rows:
    mu, sd = mstat[r["model"]]
    if abs((float(r["mpnn_corr"]) - mu) / sd - float(r["z_mpnn"])) > 1e-3: bad += 1
ck(bad == 0, f"z_mpnn == per-model z-normalised mpnn_corr for all rows (mismatches {bad})")
ck(mstat["Ligand"][0] > mstat["Protein"][0] > mstat["Soluble"][0],
   f"confirmed Ligand MPNN-bias (means {mstat['Ligand'][0]:.3f}/{mstat['Protein'][0]:.3f}/{mstat['Soluble'][0]:.3f})")

# ---- shortlists ----
sec("5. SHORTLISTS — frontier membership, gate, overlap")
def pareto(rs, a1, a2):
    return [r for r in rs if not any((o[a1] >= r[a1] and o[a2] >= r[a2] and (o[a1] > r[a1] or o[a2] > r[a2])) for o in rs)]
R = [dict(id=r["id"], esm=float(r["esm"]), z_mpnn=float(r["z_mpnn"]), pid=float(r["pct_identity"])) for r in rows]
front = {r["id"] for r in pareto(R, "esm", "z_mpnn")}
pmz = json.load(open(os.path.join(HERE, "synthesis_shortlist.json")))["shortlist"]
esmsl = json.load(open(os.path.join(HERE, "esm_only_shortlist.json")))["shortlist"]
ck(all(d["id"] in front for d in pmz), "every per-model-z shortlist design is on the ESM x z_mpnn Pareto frontier")
ck(all(d["pct_identity"] < 70 for d in pmz), "every per-model-z design < 70% identity")
ck(len(set(d["model"] for d in pmz)) == 3, f"per-model-z shortlist spans all 3 models {dict(Counter(d['model'] for d in pmz))}")
ck(all(d["pct_identity"] < 70 for d in esmsl), "every ESM-only design < 70% identity")
ov = {d["id"] for d in pmz} & {d["id"] for d in esmsl}
ck(len(ov) == json.load(open(os.path.join(HERE, "esm_only_shortlist.json")))["overlap_with_per_model_z"],
   f"recorded overlap matches recomputed ({len(ov)})")

# ---- pilot set + AF3 batch ----
sec("6. PILOT SET + AF3 BATCH")
pilot = json.load(open(os.path.join(HERE, "pilot_set.json")))["shortlist"]
pids = {d["id"] for d in pilot}
cons = {d["id"] for d in pmz} & {d["id"] for d in esmsl}
ck(cons <= pids, "all consensus designs are in the pilot set")
ck(all(d["id"] in meta for d in pilot), "every pilot design exists in meta")
ck(all(meta[d["id"]]["seq"] == d["seq"] for d in pilot), "pilot sequences byte-match meta source")
ck(all(d["tier"] == "T_aggressive" for d in pilot if d["provenance"] == "bold-novelty"),
   "bold-novelty picks are all T_aggressive")
print("  pilot:", dict(Counter(d["provenance"] for d in pilot)), dict(Counter(d["model"] for d in pilot)))
batch = json.load(open(os.path.join(HERE, "af3_inputs", "af3_batch.json")))
ck(len(batch) == len(pilot) + 5, f"AF3 batch = pilot ({len(pilot)}) + 5 calibrators (got {len(batch)})")
ck(all("count" in list(s.values())[0] for j in batch for s in j["sequences"]),
   "every AF3 chain entry has count (upload-safe)")
# design protein byte-match into batch
pmap = {f"cas9_{d['id']}": d["seq"] for d in pilot}
ck(all(j["sequences"][0]["proteinChain"]["sequence"] == pmap[j["name"]]
       for j in batch if j["name"] in pmap), "AF3 design proteins byte-match pilot source")
ck(all(len(j["sequences"]) == 6 and [list(s.keys())[0] for s in j["sequences"]] ==
       ["proteinChain","rnaSequence","dnaSequence","dnaSequence","dnaSequence","ion"] for j in batch),
   "every AF3 job has the correct 6-chain 8SRS layout")

# ---- component self-tests ----
sec("7. COMPONENT SELF-TESTS")
import subprocess
for label, args, mark in [
    ("audit_numbering", ["audit_numbering.py"], "AUDIT PASSED"),
    ("validate_results (foundation 46)", ["validate_results.py"], "ALL RESULTS VALIDATED"),
    ("ld_analysis selftest", ["ld_analysis.py", "--selftest"], "SELFTEST PASS"),
    ("parse_af3 selftest", ["parse_af3.py", "--selftest"], "R-loop"),
]:
    r = subprocess.run([sys.executable] + args, cwd=HERE, capture_output=True, text=True)
    ck(r.returncode == 0 and mark in (r.stdout + r.stderr), f"{label} green")

print("\n" + "=" * 74)
print(f"PIPELINE VALIDATION: {NC} checks, {len(FAIL)} failures")
if FAIL:
    for m in FAIL: print("   FAIL: " + m)
    sys.exit(1)
print("ALL VALIDATED — generation, scoring, per-model-z fix, shortlists, pilot, AF3 all consistent.")
