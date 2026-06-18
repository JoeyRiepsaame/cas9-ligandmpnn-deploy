#!/usr/bin/env python3
"""Parse SpRY-Cas9 AF3-Server results -> interface + active-site table (Phase 6).

Chain order (from build_af3.py): [0]Cas9 [1]sgRNA [2]DNA_C [3]DNA_D [4]DNA_c [5]MG
  Cas9:sgRNA = chain_pair_iptm[0][1]   Cas9:DNA = max [0][2..4]   R-loop = max [1][2..4]
AF3 is a FOLD/BINDING gate, not a function ranker — read it against the calibrator
ladder (WT~=dcas9 >> scram; WT >> bind_dead if AF3 sees binding loss).

Also reports per-residue pLDDT at the 7 catalytic residues + the omit-C site C80
+ a PI-domain SpRY residue (active-site / PAM-domain ordering sanity). Author
resnum -> AF3 chain-A index via wt_resnum_to_index.json (AF3 renumbers 1..N).

Usage: python parse_af3.py <unzipped_results_dir>   [--selftest]
"""
import argparse, glob, json, os, statistics as st
HERE = os.path.dirname(os.path.abspath(__file__))

r2i = {int(k): v for k, v in json.load(open(os.path.join(HERE, "wt_resnum_to_index.json"))).items()}
# author resnum -> AF3 auth_seq_id (1-based index in the design's chain A)
AF3 = {rn: i + 1 for rn, i in r2i.items()}
PROBE = {**{p: f"cat_{p}" for p in (10, 762, 839, 840, 863, 983, 986)},
         80: "exCys_C80", 1135: "SpRY_1135"}


def best_summary(job_dir):
    sums = glob.glob(os.path.join(job_dir, "*summary_confidences*.json"))
    best = None
    for s in sums:
        d = json.load(open(s))
        rk = d.get("ranking_score", d.get("iptm", 0))
        if best is None or rk > best[0]:
            model = s.replace("summary_confidences", "model").replace(".json", ".cif")
            full = s.replace("summary_confidences", "full_data")
            best = (rk, d, model if os.path.exists(model) else None,
                    full if os.path.exists(full) else None)
    return (best[1], best[2], best[3]) if best else (None, None, None)


def ca_plddt(model_path):
    """{auth_seq_id: pLDDT} for protein CA atoms (chain A), from B-factor column."""
    if not model_path or not os.path.exists(model_path):
        return {}
    out, cols, in_loop = {}, {}, False
    for line in open(model_path):
        t = line.strip()
        if t == "loop_":
            cols, in_loop = {}, True; continue
        if in_loop and t.startswith("_atom_site."):
            cols[t.split(".", 1)[1]] = len(cols); continue
        if cols and t.startswith(("ATOM", "HETATM")):
            f = t.split()
            if len(f) < len(cols):
                continue
            if f[cols["label_atom_id"]].strip('"') != "CA":
                continue
            # protein chain only (first chain); use auth_asym_id if present
            try:
                rid = int(f[cols.get("auth_seq_id", cols.get("label_seq_id"))])
                ch = f[cols["auth_asym_id"]] if "auth_asym_id" in cols else "A"
                if ch in ("A",):
                    out[rid] = float(f[cols["B_iso_or_equiv"]])
            except (ValueError, KeyError):
                pass
        elif cols and t and not t.startswith(("_", "ATOM", "HETATM")):
            in_loop = False
    return out


def mean_plddt(full_path):
    if not full_path or not os.path.exists(full_path):
        return None
    d = json.load(open(full_path)); p = d.get("atom_plddts") or d.get("plddt")
    return round(st.mean(p), 1) if p else None


def parse_dir(root):
    man = {}
    mp = os.path.join(HERE, "af3_inputs", "af3_manifest.json")
    if os.path.exists(mp):
        man = {j["name"]: j for j in json.load(open(mp))["jobs"]}
    rows = []
    for jd in sorted(glob.glob(os.path.join(root, "*"))):
        if not os.path.isdir(jd):
            continue
        name = os.path.basename(jd).lower().replace("fold_", "")
        summ, model, full = best_summary(jd)
        if summ is None:
            continue
        cpi = summ.get("chain_pair_iptm")
        def pr(i, j):
            try: return round(cpi[i][j], 2)
            except Exception: return None
        cas9_rna = pr(0, 1)
        cas9_dna = [pr(0, k) for k in (2, 3, 4) if pr(0, k) is not None]
        rloop = [pr(1, k) for k in (2, 3, 4) if pr(1, k) is not None]
        ca = ca_plddt(model)
        probe = {lab: (round(ca[AF3[rn]], 1) if AF3.get(rn) in ca else None)
                 for rn, lab in PROBE.items()}
        role = man.get("cas9_" + name, man.get(name, {})).get("role", "")
        rows.append(dict(job=name, role=role,
                         rank=round(summ.get("ranking_score", 0), 3),
                         iptm=summ.get("iptm"), ptm=summ.get("ptm"),
                         plddt=mean_plddt(full),
                         cas9_sgRNA=cas9_rna,
                         cas9_DNA=round(max(cas9_dna), 2) if cas9_dna else None,
                         rloop=round(max(rloop), 2) if rloop else None,
                         **probe))
    # calibrators first, then by R-loop desc
    rows.sort(key=lambda r: (not r["role"].startswith("calib"), -(r["rloop"] or 0)))
    hdr = f"{'job':22}{'role':16}{'iptm':>6}{'pTM':>6}{'pLDDT':>7}{'Cas9:gRNA':>10}{'Cas9:DNA':>9}{'R-loop':>8}{'cat_min':>8}"
    print(hdr); print("-" * len(hdr))
    for r in rows:
        catvals = [r[k] for k in r if k.startswith("cat_") and r[k] is not None]
        cmin = min(catvals) if catvals else None
        def f(x, w, p=2): return f"{x:>{w}.{p}f}" if isinstance(x, (int, float)) else f"{'-':>{w}}"
        print(f"{r['job']:22}{r['role']:16}{f(r['iptm'],6)}{f(r['ptm'],6)}{f(r['plddt'],7,1)}"
              f"{f(r['cas9_sgRNA'],10)}{f(r['cas9_DNA'],9)}{f(r['rloop'],8)}{f(cmin,8,1)}")
    out = os.path.join(root, "af3_interface_table.csv")
    import csv
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"\nwrote {out}")
    print("Gate reading: designs scoring in the WT/dcas9 R-loop band PASS; near scram_rec = reject.")


def selftest():
    import tempfile
    T = tempfile.mkdtemp()
    jd = os.path.join(T, "fold_cas9_d00102"); os.makedirs(jd)
    cpi = [[0.9, 0.61, 0.40, 0.55, 0.20, 0.0],
           [0.61, 0.9, 0.30, 0.72, 0.10, 0.0],
           [0.40, 0.30, 0.9, 0.25, 0.10, 0.0],
           [0.55, 0.72, 0.25, 0.9, 0.15, 0.0],
           [0.20, 0.10, 0.10, 0.15, 0.9, 0.0],
           [0.0]*6]
    json.dump({"ranking_score": 0.8, "iptm": 0.77, "ptm": 0.83, "chain_pair_iptm": cpi},
              open(os.path.join(jd, "fold_cas9_d00102_summary_confidences_0.json"), "w"))
    json.dump({"atom_plddts": [90.0, 88.0]}, open(os.path.join(jd, "fold_cas9_d00102_full_data_0.json"), "w"))
    # tiny CIF: CA atoms at AF3 index for catalytic D10 (idx = AF3[10]) and C80 site
    i10, i80 = AF3[10], AF3[80]
    cif = ("data_x\nloop_\n_atom_site.group_PDB\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n"
           "_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.B_iso_or_equiv\n"
           f"ATOM CA ASP A {i10} 95.0\nATOM CA ALA A {i80} 62.0\n")
    open(os.path.join(jd, "fold_cas9_d00102_model_0.cif"), "w").write(cif)
    print("=== parse_af3 selftest ===")
    parse_dir(T)
    print("expected: R-loop = max(cp[1][2..4]) = 0.72 ; Cas9:gRNA=0.61 ; cat_10~95 exCys_C80~62")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("results_dir", nargs="?")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest or not a.results_dir:
        selftest()
    else:
        parse_dir(a.results_dir)
