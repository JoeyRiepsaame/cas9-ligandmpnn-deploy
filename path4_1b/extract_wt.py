#!/usr/bin/env python3
"""Extract the WT SpRY-Cas9 design reference (8SRS chain A) to a FASTA + audit it.

This is the reference every downstream score uses (ESM naturalness, %identity).
Per the sequence-safety rule we generate it from the structure file via a script
and AUDIT it byte-for-byte against the structure's functional anchors, rather than
typing it. Two artifacts:

  spry_cas9_wt.fasta            modeled chain-A residues in author order (design len)
  wt_resnum_to_index.json       {author_resnum: 0-based index in the FASTA}

The audit asserts the 7 catalytic residues sit at their expected author resnums
and the 11 SpRY positions are present — if the structure parse were wrong, this
fails loudly before any sequence is used.
"""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from constraints.mapping import parse_structure_chain  # {author_resnum: one_letter}

dm = json.load(open(os.path.join(HERE, "domain_map.json")))
CATALYTIC = {c["pos"]: c["label"] for c in dm["catalytic"]}   # 10:D10, 762:E762, ...
SPRY = [s["pos"] for s in dm["spry_positions"]]

res = parse_structure_chain(os.path.join(REPO, "inputs/8SRS.cif"), chain="A")
resnums = sorted(res)                       # ascending author numbering (gaps ok)
seq = "".join(res[r] for r in resnums)      # modeled residues, author order
r2i = {r: i for i, r in enumerate(resnums)}

print(f"chain A modeled residues : {len(seq)}")
print(f"author range             : {resnums[0]}-{resnums[-1]} "
      f"(gaps = {resnums[-1]-resnums[0]+1-len(resnums)})")

# ---- AUDIT against functional anchors ----
print("\n=== catalytic residue audit (expected AA at author resnum) ===")
ok = True
for pos, label in sorted(CATALYTIC.items()):
    want = label[0]                          # 'D' from 'D10'
    got = res.get(pos, "?")
    flag = "OK" if got == want else "MISMATCH"
    ok &= (got == want)
    print(f"  {label:6} resnum {pos:>4}: structure={got}  expected={want}  {flag}")
assert ok, "CATALYTIC RESIDUE MISMATCH — structure parse or numbering is wrong, ABORT"

print("\n=== SpRY positions present ===")
missing = [p for p in SPRY if p not in res]
spry_str = ", ".join("{}:{}".format(p, res.get(p, "-")) for p in SPRY)
print(f"  {len(SPRY)-len(missing)}/{len(SPRY)} present; "
      f"structure residues at SpRY sites: {{{spry_str}}}")
assert not missing, f"SpRY positions missing from structure: {missing}"

# ---- write artifacts ----
with open(os.path.join(HERE, "spry_cas9_wt.fasta"), "w") as fh:
    fh.write(">SpRY_Cas9_8SRS_chainA_design_ref\n")
    for i in range(0, len(seq), 60):
        fh.write(seq[i:i+60] + "\n")
json.dump({str(k): v for k, v in r2i.items()},
          open(os.path.join(HERE, "wt_resnum_to_index.json"), "w"))

# round-trip check: re-read FASTA and confirm it matches
back = "".join(l.strip() for l in open(os.path.join(HERE, "spry_cas9_wt.fasta"))
               if not l.startswith(">"))
assert back == seq, "FASTA round-trip mismatch!"
print(f"\nAUDIT PASS. wrote spry_cas9_wt.fasta (len {len(seq)}) + wt_resnum_to_index.json")
print("NOTE: sequence is in the FILE; not printed here per sequence-safety rule.")
