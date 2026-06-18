#!/usr/bin/env python3
"""SpRY-Cas9 (8SRS chain A) domain map + per-domain fixation budget.

This is the FOUNDATION spec for the Path-4.1b port onto SpRY-Cas9. It encodes:
  1. SpCas9 domain boundaries (author numbering, == 8SRS chain A numbering).
  2. The project's *audited* functional anchors (NOT from memory):
       - SPRY_POSITIONS  (11) -> constraints/pipeline.py SPRY_POSITIONS
       - CATALYTIC       (7)  -> tier1_fixed.txt (A10/A762/A839/A840/A863/A983/A986)
  3. Per-domain fixation policy implementing the chosen "per-domain identity"
     strategy: PROTECT the PI/PAM-interacting domain (carries 10/11 SpRY muts +
     the PAM readout) and HNH (allosteric/RamR-risk); REDESIGN the REC lobe
     (main novelty source); keep RuvC/BH scaffold MODERATE (catalytic + fold core).
  4. The three catalytic tracks (active / nickase / dCas9) and how they are
     produced (nickase + dCas9 are POST-HOC point mutations on the active
     designs, not separate MPNN runs).

Domain boundaries are the standard SpCas9 architecture (Nishimasu 2014, Cell
156:935; Anders 2014). They are approximate at the linkers (+/-~10 aa) but the
strategy is robust to that because the functional anchors (SpRY in PI, catalytic
triads) are exact. Run `python domain_map.py` to emit domain_map.json + an audit.
"""
import json, os

# --- 1. SpCas9 domain boundaries (inclusive, author numbering = 8SRS chain A) ---
# lobe membership noted; RuvC is a split domain (I/II/III) interleaved with HNH.
DOMAINS = [
    # name            start  end    lobe     policy
    ("RuvC-I",            1,    59,  "NUC",   "MODERATE"),
    ("BridgeHelix",     60,    93,  "NUC",   "MODERATE"),  # Arg-rich; A61 (SpRY) here
    ("REC1",            94,   179,  "REC",   "REDESIGN"),
    ("REC2",           180,   307,  "REC",   "REDESIGN"),
    ("REC3",           308,   717,  "REC",   "REDESIGN"),
    ("RuvC-II",        718,   774,  "NUC",   "MODERATE"),
    ("HNH",            775,   908,  "NUC",   "PROTECT"),   # allosteric/catalytic
    ("RuvC-III",       909,  1098,  "NUC",   "MODERATE"),
    ("PI",            1099,  1368,  "NUC",   "PROTECT"),   # PAM-interacting; 10/11 SpRY
]

# --- 2. AUDITED functional anchors (sourced from the repo, not memory) ---
# constraints/pipeline.py : SPRY_POSITIONS
SPRY_POSITIONS = [61, 1111, 1135, 1136, 1218, 1219, 1317, 1322, 1333, 1335, 1337]
# tier1_fixed.txt catalytic set (RuvC: D10,E762,H983,D986 ; HNH: D839,H840,N863)
CATALYTIC = {
    10:  ("D10",  "RuvC-I"),
    762: ("E762", "RuvC-II"),
    839: ("D839", "HNH"),
    840: ("H840", "HNH"),
    863: ("N863", "HNH"),
    983: ("H983", "RuvC-III"),
    986: ("D986", "RuvC-III"),
}

# --- 3. Per-domain fixation policy (the "per-domain identity" strategy) ---
# target_identity = rough per-domain identity-to-WT we steer toward via tier choice.
# conservation_floor = the conservation threshold whose fixed set seeds this domain.
POLICY = {
    "PROTECT":  dict(target_identity="high (>=70%)", conservation_floor=50,
                     note="high fixation: keep function-defining surface (PAM / HNH catalysis+docking)"),
    "MODERATE": dict(target_identity="mid (~60-70%)", conservation_floor=70,
                     note="catalytic residues + fold core fixed; scaffold redesignable"),
    "REDESIGN": dict(target_identity="low (<70%, novelty)", conservation_floor=90,
                     note="main novelty source: only the most-conserved core fixed"),
}

# --- 4. Catalytic tracks (nickase + dCas9 are post-hoc point mutations) ---
TRACKS = {
    "active": dict(mutations={}, note="all 7 catalytic WT; functional dsDNA cleavage"),
    "nickase_nCas9": dict(mutations={840: "A"},
                          note="H840A: HNH-dead nickase (canonical for base/prime editing)"),
    "dCas9": dict(mutations={10: "A", 840: "A"},
                  note="D10A+H840A: catalytically dead (CRISPRi/anchoring)"),
}
# IMPORTANT: design ONCE with all 7 catalytic fixed-WT (active track). Derive
# nickase/dCas9 by installing the point mutation(s) afterward -> ~0 extra MPNN cost
# and preserves the catalytic geometry MPNN conditions on.


def domain_of(pos):
    for name, lo, hi, lobe, pol in DOMAINS:
        if lo <= pos <= hi:
            return name, lobe, pol
    return None, None, None


def build():
    domains = [dict(name=n, start=lo, end=hi, length=hi - lo + 1, lobe=lobe,
                    policy=pol, **POLICY[pol]) for (n, lo, hi, lobe, pol) in DOMAINS]
    spry = [dict(pos=p, domain=domain_of(p)[0]) for p in SPRY_POSITIONS]
    catalytic = [dict(pos=p, label=lbl, domain_expected=dom, domain_computed=domain_of(p)[0])
                 for p, (lbl, dom) in sorted(CATALYTIC.items())]
    return dict(
        target="SpRY-Cas9 / 8SRS chain A",
        numbering="author (8SRS chain A, 3-1366, 23 internal gaps)",
        full_length_ref=1368,
        strategy="per-domain identity: PROTECT PI+HNH, REDESIGN REC, MODERATE RuvC/BH",
        domains=domains,
        spry_positions=spry,
        catalytic=catalytic,
        tracks=TRACKS,
        mandatory_locks=sorted(set(SPRY_POSITIONS) | set(CATALYTIC)),
    )


def audit(m):
    print("=== DOMAIN MAP AUDIT ===")
    # a) every SpRY position lands in PI or BridgeHelix
    bad = [s for s in m["spry_positions"] if s["domain"] not in ("PI", "BridgeHelix")]
    pi = sum(s["domain"] == "PI" for s in m["spry_positions"])
    bh = sum(s["domain"] == "BridgeHelix" for s in m["spry_positions"])
    print(f"SpRY: {len(m['spry_positions'])} positions -> PI={pi}, BridgeHelix={bh}, other={len(bad)}")
    assert not bad, f"SpRY positions outside PI/BH: {bad}"
    # b) catalytic domain expected == computed
    mism = [c for c in m["catalytic"] if c["domain_expected"] != c["domain_computed"]]
    print(f"Catalytic: {len(m['catalytic'])} residues, domain mismatches={len(mism)}")
    assert not mism, f"catalytic domain mismatch: {mism}"
    # c) domain coverage / no overlaps / lengths
    spans = sorted((d["start"], d["end"], d["name"]) for d in m["domains"])
    for (s1, e1, n1), (s2, e2, n2) in zip(spans, spans[1:]):
        assert e1 < s2, f"overlap/gap between {n1} and {n2}"
    redesign = sum(d["length"] for d in m["domains"] if d["policy"] == "REDESIGN")
    protect  = sum(d["length"] for d in m["domains"] if d["policy"] == "PROTECT")
    moderate = sum(d["length"] for d in m["domains"] if d["policy"] == "MODERATE")
    tot = redesign + protect + moderate
    print(f"Residue budget (of {tot} mapped): REDESIGN(REC)={redesign} "
          f"({100*redesign/tot:.0f}%) | PROTECT(PI+HNH)={protect} ({100*protect/tot:.0f}%) "
          f"| MODERATE(RuvC/BH)={moderate} ({100*moderate/tot:.0f}%)")
    print("ALL CHECKS PASS\n")


if __name__ == "__main__":
    m = build()
    audit(m)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "domain_map.json")
    json.dump(m, open(out, "w"), indent=2)
    print(f"wrote {out}")
    print("\nPer-domain plan:")
    for d in m["domains"]:
        print(f"  {d['name']:12} {d['start']:>4}-{d['end']:<4} {d['lobe']:4} "
              f"{d['policy']:9} -> identity {d['target_identity']}")
