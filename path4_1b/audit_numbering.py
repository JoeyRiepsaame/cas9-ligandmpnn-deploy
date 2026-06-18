#!/usr/bin/env python3
"""THOROUGH numbering + pipeline audit for the SpRY-Cas9 Path-4.1b campaign.

The whole campaign rides on one assumption: the author residue numbers we fix
(catalytic, SpRY, conservation, DCA) really are canonical SpCas9 positions, and
the gap-safe MSA->PDB zip in constraints/mapping.py did not silently slip after
an unmodeled loop. This script proves it three independent ways:

  A. EXTERNAL alignment  : align the 8SRS chain-A WT (1341 modeled) to canonical
                           full-length SpCas9 (UniProt Q99ZW2, 1368) and assert
                           every modeled residue's author resnum == its aligned
                           canonical position. (ground-truth numbering)
  B. INTERNAL consistency: compare the MSA's WT query (1341 ungapped) against the
                           extracted 8SRS WT index-by-index; the differences must
                           be EXACTLY the 11 SpRY positions (and nothing else).
  C. PIPELINE reproduce  : re-run build_column_to_pdb on the real MSA (0 unexpected
                           mismatches) and confirm conservation/DCA author resnums
                           all map onto real modeled residues.

Plus: characterize the SpRY substitutions (canonical->structure) and confirm the
catalytic residues are untouched by the engineering. Exit non-zero on any failure.
"""
import gzip, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from constraints.mapping import parse_structure_chain, build_column_to_pdb
from Bio import Align

FAIL = []
def check(cond, msg):
    print(("  [PASS] " if cond else "  [FAIL] ") + msg)
    if not cond:
        FAIL.append(msg)

def read_fasta(path, opener=open):
    seq = []
    with opener(path, "rt") as f:
        for ln in f:
            if not ln.startswith(">"):
                seq.append(ln.strip())
    return "".join(seq)

# ---------- load everything ----------
struct = parse_structure_chain(os.path.join(REPO, "inputs/8SRS.cif"), chain="A")
pdb_sorted = sorted(struct)                       # ascending author resnums (1341)
wt8srs = "".join(struct[r] for r in pdb_sorted)   # modeled residues, author order
canon = read_fasta(os.path.join(HERE, "ref/Q99ZW2_SpCas9_canonical.fasta"))
dm = json.load(open(os.path.join(HERE, "domain_map.json")))
SPRY = sorted(s["pos"] for s in dm["spry_positions"])
CATALYTIC = {c["pos"]: c["label"][0] for c in dm["catalytic"]}

# MSA WT query (row 0), ungapped
with gzip.open(os.path.join(REPO, "constraints/data/cas9_msa.fasta.gz"), "rt") as f:
    f.readline()                                  # header
    msa_wt_aligned = ""
    for ln in f:
        if ln.startswith(">"):
            break
        msa_wt_aligned += ln.strip()
msa_wt = msa_wt_aligned.replace("-", "")

print(f"8SRS chain A modeled : {len(wt8srs)}  (author {pdb_sorted[0]}-{pdb_sorted[-1]})")
print(f"canonical SpCas9     : {len(canon)}")
print(f"MSA WT query ungapped: {len(msa_wt)}\n")

# =================== A. EXTERNAL alignment vs canonical SpCas9 ===================
print("A. EXTERNAL alignment audit (author resnum == canonical position)")
al = Align.PairwiseAligner()
al.mode = "global"; al.open_gap_score = -11; al.extend_gap_score = -1
al.match_score = 5; al.mismatch_score = -4
aln = al.align(canon, wt8srs)[0]                  # target=canonical, query=8SRS
# walk aligned blocks: canonical[cs:ce] <-> query[qs:qe]
num_mismatch = 0; subs = []; mapped = 0
for (cs, ce), (qs, qe) in zip(aln.aligned[0], aln.aligned[1]):
    for k in range(qe - qs):
        qi = qs + k                               # index into wt8srs
        canon_pos = cs + k + 1                    # 1-based canonical position
        author = pdb_sorted[qi]                   # author resnum of this modeled residue
        mapped += 1
        if author != canon_pos:
            num_mismatch += 1
            if len(subs) < 10:
                subs.append(f"author{author}!=canon{canon_pos}")
        if wt8srs[qi] != canon[canon_pos - 1]:
            subs_aa = (canon_pos, canon[canon_pos - 1], wt8srs[qi])
check(mapped == len(wt8srs), f"all {len(wt8srs)} modeled residues placed by alignment (got {mapped})")
check(num_mismatch == 0, f"author resnum == canonical position for ALL modeled residues "
                         f"(offsets={num_mismatch}{'; '+', '.join(subs) if subs else ''})")
# unmodeled residues decompose as: internal gaps (within author range) + terminal residues
lo, hi = pdb_sorted[0], pdb_sorted[-1]
internal_gaps = (hi - lo + 1) - len(wt8srs)        # missing within [lo,hi]
terminal = (lo - 1) + (len(canon) - hi)            # canon residues outside [lo,hi]
total_unmodeled = len(canon) - len(wt8srs)
check(internal_gaps == 23, f"internal gaps within author {lo}-{hi} == 23 (got {internal_gaps})")
check(terminal == 4, f"unmodeled terminal residues == 4 "
                     f"(N-term 1..{lo-1}, C-term {hi+1}..{len(canon)}; got {terminal})")
check(internal_gaps + terminal == total_unmodeled,
      f"23 internal + 4 terminal == {total_unmodeled} canonical-minus-modeled")

# =================== B. INTERNAL: MSA WT query vs 8SRS WT ===================
print("\nB. INTERNAL consistency (MSA WT query vs structure; diffs == SpRY only)")
check(len(msa_wt) == len(wt8srs), f"MSA WT query and structure same length ({len(msa_wt)}=={len(wt8srs)})")
diffs = [pdb_sorted[i] for i in range(min(len(msa_wt), len(wt8srs))) if msa_wt[i] != wt8srs[i]]
check(set(diffs) == set(SPRY),
      f"index-by-index diffs are EXACTLY the 11 SpRY positions "
      f"(got {len(diffs)}: {diffs[:15]})")

# =================== C. PIPELINE reproduce (mapping + conservation/DCA) ===================
print("\nC. PIPELINE reproduction (mapping 0-mismatch; conservation/DCA mapped)")
col_to_pdb, report = build_column_to_pdb(msa_wt_aligned, struct, spry_positions=SPRY, strict=False)
check(not report["unexpected_mismatches"],
      f"build_column_to_pdb: 0 non-SpRY mismatches (got {len(report['unexpected_mismatches'])})")
check(report["columns_mapped"] == len(wt8srs),
      f"columns_mapped == modeled residues ({report['columns_mapped']})")
# every conservation-tier + DCA author resnum must be a real modeled residue
con = json.load(open(os.path.join(REPO, "outputs/constraints/8srs_constraints.json")))
modeled = set(pdb_sorted)
allcons = set().union(*[set(con["conservation_tiers"][k]) for k in con["conservation_tiers"]])
dca = set(con["dca_top"])
check(allcons <= modeled, f"all conservation-tier positions are modeled residues "
                          f"(stray={sorted(allcons-modeled)[:5]})")
check(dca <= modeled, f"all {len(dca)} DCA positions are modeled residues "
                      f"(stray={sorted(dca-modeled)[:5]})")
# spot-check the mapping at a catalytic anchor: the column mapping to pdb 10 must have WT query D
col_of_10 = [c for c, p in col_to_pdb.items() if p == 10]
check(len(col_of_10) == 1 and msa_wt_aligned[col_of_10[0]] == "D",
      f"MSA column -> pdb resnum 10 carries WT query 'D' (catalytic D10 anchor)")

# =================== D. SpRY substitution + catalytic characterization ===================
print("\nD. SpRY substitutions (canonical -> 8SRS) + catalytic integrity")
known_spry = {61:("A","R"),1111:("L","R"),1135:("D","L"),1136:("S","W"),1218:("G","K"),
              1219:("E","Q"),1317:("N","R"),1322:("A","R"),1333:("R","P"),1335:("R","Q"),1337:("T","R")}
ok_spry = True
for p in SPRY:
    cano = canon[p-1]; got = struct[p]; want_from, want_to = known_spry[p]
    good = (cano == want_from and got == want_to)
    ok_spry &= good
    print(f"    {p:>4}: canonical {cano} -> structure {got}   "
          f"(expect {want_from}->{want_to}) {'OK' if good else 'MISMATCH'}")
check(ok_spry, "all 11 SpRY substitutions match the published SpRY mutation set")
ok_cat = all(canon[p-1] == aa and struct[p] == aa for p, aa in CATALYTIC.items())
check(ok_cat, "all 7 catalytic residues identical in canonical AND structure (engineering left them WT)")

# ---------- verdict ----------
print("\n" + "=" * 64)
if FAIL:
    print(f"AUDIT FAILED — {len(FAIL)} check(s):")
    for m in FAIL:
        print("   - " + m)
    sys.exit(1)
print("AUDIT PASSED — numbering is canonical SpCas9 end-to-end; pipeline consistent.")
