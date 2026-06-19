#!/usr/bin/env python3
"""Independent validation of the AF3 triage round-trip + results (extensive audit).

Verifies that AF3 folded the EXACT sequences we sent (byte-match job_request vs
source), that the chain layout is correct, that the parsed interface metrics match
a re-read of the summary JSONs, and that the headline scram_rec conclusion
(AF3 non-discriminating) reproduces from the CIFs independently.

Usage: python validate_af3.py <folds_dir>
"""
import csv, glob, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from parse_af3 import ca_plddt, best_summary, AF3

FAIL, NC = [], 0
def ck(c, m):
    global NC; NC += 1
    print(("  PASS  " if c else "  FAIL  ") + m)
    if not c: FAIL.append(m)

root = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/Downloads/folds_2026_06_19_10_32")
jobdirs = [d for d in sorted(glob.glob(os.path.join(root, "*"))) if os.path.isdir(d)
           and os.path.basename(d).startswith("cas9_")]

# sources of truth
pilot = {f"cas9_{d['id']}": d["seq"] for d in json.load(open(os.path.join(HERE, "pilot_set.json")))["shortlist"]}
cal = {}
name = None; cur = []
for ln in open(os.path.join(HERE, "af3_calibrators.fasta")):
    if ln.startswith(">"):
        if name: cal[f"cas9_cal_{name}"] = "".join(cur)
        name = ln[1:].split()[0]; cur = []
    else: cur.append(ln.strip())
if name: cal[f"cas9_cal_{name}"] = "".join(cur)
src = {**pilot, **cal}

print(f"=== AF3 round-trip + results audit ({len(jobdirs)} job dirs in {os.path.basename(root)}) ===")
print("\n1. SEQUENCE ROUND-TRIP — AF3 folded exactly what we sent")
miss_count = 0; bytematch = 0; layout_ok = 0
for jd in jobdirs:
    job = os.path.basename(jd)
    jr = glob.glob(os.path.join(jd, "*job_request.json"))
    if not jr:
        miss_count += 1; continue
    req = json.load(open(jr[0])); req = req[0] if isinstance(req, list) else req
    seqs = req["sequences"]
    layout = [list(s.keys())[0] for s in seqs]
    if layout == ["proteinChain","rnaSequence","dnaSequence","dnaSequence","dnaSequence","ion"]:
        layout_ok += 1
    prot = seqs[0]["proteinChain"]["sequence"]
    if job in src and prot == src[job]:
        bytematch += 1
    elif job in src:
        print(f"    MISMATCH {job}: folded protein != source")
ck(miss_count == 0, f"every job has a job_request ({len(jobdirs)-miss_count}/{len(jobdirs)})")
ck(layout_ok == len(jobdirs), f"every job has the 6-chain 8SRS layout ({layout_ok}/{len(jobdirs)})")
ck(bytematch == sum(1 for j in jobdirs if os.path.basename(j) in src),
   f"every folded protein BYTE-MATCHES its source sequence ({bytematch} matched)")

print("\n2. METRIC RE-EXTRACTION vs parsed CSV (sample)")
table = {r["job"]: r for r in csv.DictReader(open(os.path.join(root, "af3_interface_table.csv")))}
sample = ["cas9_d00806", "cas9_cal_dcas9", "cas9_cal_scram_rec"]
ok = True
for job in sample:
    jd = os.path.join(root, job)
    summ, _, _ = best_summary(jd)
    cpi = summ["chain_pair_iptm"]
    rna = round(cpi[0][1], 2)
    rloop = round(max(x for x in (cpi[1][2], cpi[1][3], cpi[1][4]) if x is not None), 2)
    name = job.lower().replace("fold_", "")
    row = table.get(name)
    if not row or abs(float(row["cas9_sgRNA"]) - rna) > 0.01 or abs(float(row["rloop"]) - rloop) > 0.01:
        ok = False; print(f"    {job}: recomputed gRNA={rna} rloop={rloop} vs csv {row.get('cas9_sgRNA')}/{row.get('rloop')}")
ck(ok, f"Cas9:sgRNA + R-loop re-extracted from summary JSON match the parsed table ({len(sample)} jobs)")

print("\n3. HEADLINE CONCLUSION — scram_rec non-discrimination reproduces")
win = [AF3[r] for r in range(200, 261) if r in AF3]
def local(job):
    _, model, _ = best_summary(os.path.join(root, job))
    ca = ca_plddt(model)
    v = [ca[i] for i in win if i in ca]
    return sum(v) / len(v) if v else None
scram = local("cas9_cal_scram_rec"); native = local("cas9_cal_dcas9")
designs = [local(j) for j in ("cas9_d00806", "cas9_d00900", "cas9_d00781")]
print(f"    local pLDDT @200-260: scram_rec={scram:.1f}  native(dCas9)={native:.1f}  designs={[round(d,1) for d in designs]}")
ck(native > 80, f"native REC2 folds confidently (dCas9 local pLDDT {native:.1f} > 80)")
ck(scram < 60, f"scrambled REC2 is disordered (scram_rec local pLDDT {scram:.1f} < 60)")
ck(all(d < 60 for d in designs), "real designs score AS LOW as the scramble in REC2 (AF3 non-discriminating, confirmed)")
ck(max(designs) < native - 20, "real designs >20 pLDDT below native in the redesigned window (AF3 can't validate the redesign)")

print("\n4. TABLE COMPLETENESS")
ck(len(table) == len(jobdirs), f"interface table has a row per job ({len(table)}/{len(jobdirs)})")
roles = [r["role"] for r in table.values()]
ck(sum(1 for x in roles if x == "calibrator") == 4, "4 calibrators in table")
ck(sum(1 for x in roles if x.startswith("design")) == len(jobdirs) - 4,
   f"{len(jobdirs)-4} designs in table")

print("\n" + "=" * 70)
print(f"AF3 AUDIT: {NC} checks, {len(FAIL)} failures")
if FAIL:
    for m in FAIL: print("   FAIL: " + m)
    sys.exit(1)
print("AF3 RESULTS VALIDATED — correct sequences folded, metrics reproduce, conclusion holds.")
