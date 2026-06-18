#!/usr/bin/env python3
"""SpRY-Cas9 AF3 calibration panel (gap C4) — built from the WT, file-based + audited.

RT lesson: AF3 R-loop iPTM measures FOLD + nucleic-acid BINDING, NOT catalysis — it
could not rank alive-vs-dead RTs. So a naive "dead = dCas9" is a BAD AF3 negative,
because dCas9 folds and binds identically to WT (only catalysis is removed). To
make the AF3 gate interpretable we need a CALIBRATION LADDER that isolates what AF3
can and cannot see:

  role                seq                         expected AF3 (R-loop iPTM / pLDDT)
  ------------------  --------------------------  ----------------------------------
  positive (native)   WT SpRY-Cas9                HIGH  (folds + binds + cuts)
  catalysis-blind     dCas9  (D10A + H840A)       ~= WT  -> proves AF3 ignores catalysis
  nickase ref         nCas9  (H840A)              ~= WT
  binding-deficient   bridge-helix R->E cluster   LOWER R-loop IF AF3 sees binding loss
  fold-negative       scrambled REC window        LOW   -> proves AF3 penalises nonsense

Reading the ladder after folding:
  * WT ~= dCas9 >> scram        -> AF3 sees fold/binding but is blind to catalysis
                                   (RT-consistent) => use AF3 ONLY as a fold/binding gate.
  * WT >> bind_dead             -> AF3 DOES register R-loop binding disruption (a real signal).
  * if WT ~= scram              -> AF3 cannot even tell nonsense in REC apart from real
                                   => AF3 is useless as a gate for our REC-heavy redesign.
The WT/dCas9 band sets the PASS threshold; designs scoring near scram are rejected.

All calibrators are derived from spry_cas9_wt.fasta by indexed edits (no transcription)
and audited (diff == intended positions only). scram_rec + bind_dead are INTENTIONALLY
non-functional COMPUTATIONAL CONTROLS — never send them to synthesis.
"""
import json, os, random
HERE = os.path.dirname(os.path.abspath(__file__))
SPRY = {61, 1111, 1135, 1136, 1218, 1219, 1317, 1322, 1333, 1335, 1337}

wt = "".join(l.strip() for l in open(os.path.join(HERE, "spry_cas9_wt.fasta")) if not l.startswith(">"))
r2i = {int(k): v for k, v in json.load(open(os.path.join(HERE, "wt_resnum_to_index.json"))).items()}

def mutate(seq, edits):
    """edits: {author_resnum: new_aa}. Returns new sequence; asserts resnum present."""
    s = list(seq)
    for rn, aa in edits.items():
        assert rn in r2i, f"resnum {rn} not modeled"
        s[r2i[rn]] = aa
    return "".join(s)

# bridge-helix arginine cluster (60-93) — data-driven; EXCLUDE the SpRY R61 (A61R)
# so the binding-dead control does not perturb the relaxed-PAM mutation set.
bh_args = [rn for rn in range(60, 94)
           if rn in r2i and wt[r2i[rn]] == "R" and rn not in SPRY][:5]
bind_dead_edits = {rn: "E" for rn in bh_args}   # charge-reverse the heteroduplex-contacting Args

# scramble a REC2 window (REDESIGN domain, no catalytic/SpRY/contacts inside): shuffle, fixed seed
REC_WIN = (200, 260)
def scramble_rec(seq):
    i0, i1 = r2i[REC_WIN[0]], r2i[REC_WIN[1]]
    block = list(seq[i0:i1 + 1])
    random.Random(42).shuffle(block)
    return seq[:i0] + "".join(block) + seq[i1 + 1:]

calibrators = {
    "wt_spry":      dict(seq=wt, role="positive/native",
                         note="active SpRY-Cas9 (folds+binds+cuts) — positive anchor"),
    "dcas9":        dict(seq=mutate(wt, {10: "A", 840: "A"}), role="catalysis-blind",
                         note="D10A+H840A: dead nuclease, fold/bind intact -> tests AF3 catalysis-blindness"),
    "ncas9_h840a":  dict(seq=mutate(wt, {840: "A"}), role="nickase ref",
                         note="H840A nickase reference"),
    "bind_dead":    dict(seq=mutate(wt, bind_dead_edits), role="binding-deficient (predicted)",
                         note=f"bridge-helix Arg->Glu at {bh_args}: predicted R-loop binding loss"),
    "scram_rec":    dict(seq=scramble_rec(wt), role="fold-negative (control)",
                         note=f"REC2 window {REC_WIN} shuffled (seed42), composition preserved -> should misfold"),
}

# ---------------- AUDIT (diff vs WT == intended positions only) ----------------
print("=== AF3 calibrator audit (diff vs WT) ===")
expected = {
    "wt_spry": set(),
    "dcas9": {10, 840},
    "ncas9_h840a": {840},
    "bind_dead": set(bh_args),
    "scram_rec": None,   # window-bounded, checked separately
}
i2r = {v: k for k, v in r2i.items()}
ok = True
for name, d in calibrators.items():
    s = d["seq"]
    assert len(s) == len(wt), f"{name}: length changed!"
    diff = {i2r[i] for i in range(len(wt)) if s[i] != wt[i]}
    if name == "scram_rec":
        inside = all(REC_WIN[0] <= rn <= REC_WIN[1] for rn in diff)
        good = inside and len(diff) > 0
        print(f"  {name:12} role={d['role']:26} changed {len(diff)} pos, all within REC{REC_WIN}: {inside}")
    else:
        exp = expected[name]
        good = (diff == exp)
        print(f"  {name:12} role={d['role']:26} diff={sorted(diff)} expected={sorted(exp)} {'OK' if good else 'MISMATCH'}")
    # safety: NO calibrator (control or reference) may perturb a SpRY position
    touched_spry = {rn for rn in diff if rn in SPRY}
    assert not touched_spry, f"{name}: altered SpRY positions {touched_spry}!"
    ok &= good
assert ok, "calibrator diff audit FAILED"

# ---------------- write artifacts ----------------
with open(os.path.join(HERE, "af3_calibrators.fasta"), "w") as fh:
    for name, d in calibrators.items():
        fh.write(f">{name} role={d['role'].replace(' ','_')}\n")
        for i in range(0, len(d["seq"]), 60):
            fh.write(d["seq"][i:i+60] + "\n")
manifest = {name: {"role": d["role"], "note": d["note"], "length": len(d["seq"]),
                   "synthesize_safe": d["role"] in ("positive/native", "nickase ref")}
            for name, d in calibrators.items()}
manifest["_meta"] = {"bridge_helix_args_mutated": bh_args, "rec_scramble_window": list(REC_WIN),
                     "intended_nonfunctional": ["bind_dead", "scram_rec"]}
json.dump(manifest, open(os.path.join(HERE, "af3_calibrators_manifest.json"), "w"), indent=2)
print(f"\nAUDIT PASS. wrote af3_calibrators.fasta ({len(calibrators)} seqs) + manifest.")
print("Include these in the Phase-6 AF3 batch alongside the design shortlist; they")
print("calibrate the R-loop iPTM PASS band. scram_rec/bind_dead are controls — NOT for synthesis.")
