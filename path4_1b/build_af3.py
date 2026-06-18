#!/usr/bin/env python3
"""Build the SpRY-Cas9 AF3-Server triage batch (Phase 6).

One AF3 job per protein = Cas9 (design or calibrator) + the 8SRS nucleic-acid
context (sgRNA + 3 DNA strands) + catalytic Mg. AF3 is a FOLD/BINDING sanity gate
only (it cannot rank catalysis or test PAM relaxation) — the calibrator ladder
(af3_calibrators) sets the interpretable PASS band.

Chain order (FIXED — parse_af3.py depends on it; chain_pair_iptm is in this order):
  [0] Cas9 protein   [1] sgRNA(chain B,98)   [2] DNA C(13)   [3] DNA D(19)   [4] DNA c(19)
  => Cas9:sgRNA = cp[0][1] · Cas9:DNA = max cp[0][2..4] · R-loop = max cp[1][2..4]

Sequences: protein from synthesis_shortlist.json + af3_calibrators.fasta (verified
files); nucleic acids parsed from 8SRS.cif. Every protein seq byte-checked; every
chain entry carries count=1 (AF3 rejects missing count). No hand transcription.

Usage: python build_af3.py [--n-designs 17 --seeds 1] -> af3_inputs/af3_batch.json + manifest
"""
import argparse, json, os
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

DNA3 = {"DA": "A", "DT": "T", "DG": "G", "DC": "C"}
RNA = {"A", "U", "G", "C"}

def parse_na(cif, chain):
    """nucleotide sequence (author order) for a nucleic-acid chain; waters skipped."""
    seq, seen = [], set()
    in_atom, cols, idx = False, [], {}
    for line in open(cif):
        if line.startswith("_atom_site."):
            cols.append(line.strip().split(".")[1]); in_atom = True; continue
        if in_atom and (line.startswith("ATOM") or line.startswith("HETATM")):
            if not idx: idx = {c: i for i, c in enumerate(cols)}
            p = line.split()
            if p[idx["auth_asym_id"]] != chain:
                continue
            comp = p[idx["auth_comp_id"]]; num = p[idx["auth_seq_id"]]
            if (chain, num) in seen:
                continue
            if comp in DNA3:
                seen.add((chain, num)); seq.append(DNA3[comp])
            elif comp in RNA:
                seen.add((chain, num)); seq.append(comp)
        elif in_atom and line.startswith("#"):
            in_atom = False
    return "".join(seq)


def load_proteins(n_designs, shortlist_file="synthesis_shortlist.json"):
    """ordered list of (job_name, role, sequence): designs from the shortlist + calibrators."""
    out = []
    sl = json.load(open(os.path.join(HERE, shortlist_file)))["shortlist"][:n_designs]
    for d in sl:
        prov = d.get("provenance", d["tier"])
        out.append((f"cas9_{d['id']}", f"design/{prov}", d["seq"]))
    # calibrators (wt/dcas9/ncas9/bind_dead/scram_rec)
    name = None; cur = []
    for ln in open(os.path.join(HERE, "af3_calibrators.fasta")):
        if ln.startswith(">"):
            if name: out.append((f"cas9_cal_{name}", "calibrator", "".join(cur)))
            name = ln[1:].split()[0]; cur = []
        else:
            cur.append(ln.strip())
    if name: out.append((f"cas9_cal_{name}", "calibrator", "".join(cur)))
    return out


def main(a):
    cif = os.path.join(REPO, "inputs/8SRS.cif")
    sgRNA = parse_na(cif, "B")
    dna = {c: parse_na(cif, c) for c in ("C", "D", "c")}
    print(f"context: sgRNA(B)={len(sgRNA)}nt  DNA C={len(dna['C'])} D={len(dna['D'])} c={len(dna['c'])}")
    assert set(sgRNA) <= RNA and len(sgRNA) == 98, "sgRNA parse off"

    proteins = load_proteins(a.n_designs, a.shortlist)
    # bio-safety: every protein 1341 aa (designs/calibrators) and present
    for nm, role, seq in proteins:
        assert len(seq) == 1341, f"{nm}: length {len(seq)} != 1341"

    jobs = []
    for nm, role, seq in proteins:
        jobs.append({
            "name": nm,
            "modelSeeds": list(range(1, a.seeds + 1)),
            "sequences": [
                {"proteinChain": {"sequence": seq, "count": 1}},
                {"rnaSequence": {"sequence": sgRNA, "count": 1}},
                {"dnaSequence": {"sequence": dna["C"], "count": 1}},
                {"dnaSequence": {"sequence": dna["D"], "count": 1}},
                {"dnaSequence": {"sequence": dna["c"], "count": 1}},
                {"ion": {"ion": "MG", "count": 1}},
            ],
        })
    # audit: count on every chain entry
    miss = sum(1 for j in jobs for s in j["sequences"]
               if "count" not in list(s.values())[0])
    assert miss == 0, f"{miss} chain entries missing count (AF3 upload would fail)"

    os.makedirs(os.path.join(HERE, "af3_inputs"), exist_ok=True)
    out = os.path.join(HERE, "af3_inputs", "af3_batch.json")
    json.dump(jobs, open(out, "w"), indent=2)
    manifest = [{"name": nm, "role": role, "protein_len": len(seq)} for nm, role, seq in proteins]
    json.dump({"chain_order": ["Cas9", "sgRNA", "DNA_C", "DNA_D", "DNA_c", "MG"],
               "rloop_iptm": "max chain_pair_iptm[1][2..4]",
               "cas9_sgRNA_iptm": "chain_pair_iptm[0][1]",
               "jobs": manifest},
              open(os.path.join(HERE, "af3_inputs", "af3_manifest.json"), "w"), indent=2)
    print(f"wrote {out}: {len(jobs)} jobs "
          f"({sum(1 for _,r,_ in proteins if r.startswith('design'))} designs + "
          f"{sum(1 for _,r,_ in proteins if r=='calibrator')} calibrators), {a.seeds} seed(s) each")
    print("Upload af3_batch.json to alphafoldserver.com; then parse_af3.py on the unzipped results.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-designs", type=int, default=18, dest="n_designs")
    ap.add_argument("--shortlist", default="synthesis_shortlist.json")
    ap.add_argument("--seeds", type=int, default=1)
    main(ap.parse_args())
