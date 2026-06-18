#!/usr/bin/env python3
"""
MPNN Blind Spot Detector — detect native residues forcibly mutated by
--omit_AAs / --bias_AAs and the artificial co-variation patterns they create.

Universal tool: works on any MPNN checkpoint JSON, not just MMLV-RT.

Usage:
    python mpnn_blind_spot_detector.py checkpoint.json --wt SEQUENCE [OPTIONS]

    # With omit_AAs constraint:
    python mpnn_blind_spot_detector.py checkpoint.json --wt SEQUENCE --omit_AAs C

    # With fixed positions (0-indexed):
    python mpnn_blind_spot_detector.py checkpoint.json --wt SEQUENCE --fixed 33,50,62,...

    # Full pipeline mode (reads fixation scheme JSON):
    python mpnn_blind_spot_detector.py checkpoint.json --wt SEQUENCE --scheme scheme.json --path path3_5b_noCys
"""

import argparse
import json
import sys
import math
import numpy as np
from collections import Counter


def parse_checkpoint(path):
    """Load MPNN checkpoint and extract design sequences."""
    with open(path) as f:
        data = json.load(f)

    if "all_designs" in data:
        designs = data["all_designs"]
        seqs = [d["sequence"] for d in designs]
    elif isinstance(data, list):
        seqs = [d["sequence"] if isinstance(d, dict) else d for d in data]
    else:
        raise ValueError(f"Cannot parse checkpoint format. Keys: {list(data.keys())}")

    return seqs


def identify_fixed_free(seqs, fixed_positions=None):
    """Identify fixed vs free positions from the design ensemble."""
    seq_len = len(seqs[0])

    if fixed_positions is not None:
        fixed = set(fixed_positions)
        free = [i for i in range(seq_len) if i not in fixed]
        return sorted(fixed), sorted(free)

    fixed = []
    free = []
    for i in range(seq_len):
        aas = set(s[i] for s in seqs)
        if len(aas) == 1:
            fixed.append(i)
        else:
            free.append(i)
    return fixed, free


def detect_omit_victims(seqs, wt_seq, free_positions, omit_aas):
    """
    Detect native residues at free positions that are forcibly mutated
    because they match --omit_AAs.

    Returns list of dicts with position, wt_aa, and design distributions.
    """
    if not omit_aas:
        return []

    omit_set = set(omit_aas.upper())
    victims = []

    for pos in free_positions:
        wt_aa = wt_seq[pos]
        if wt_aa.upper() not in omit_set:
            continue

        design_aas = [s[pos] for s in seqs]
        counter = Counter(design_aas)
        n_wt = counter.get(wt_aa, 0)
        n_total = len(seqs)

        victims.append({
            "position_0idx": pos,
            "position_1idx": pos + 1,
            "wt_aa": wt_aa,
            "omitted": True,
            "n_designs_with_wt": n_wt,
            "n_total": n_total,
            "wt_fraction": n_wt / n_total,
            "distribution": dict(counter.most_common()),
            "forced_departure": n_wt == 0,
        })

    return victims


def detect_bias_victims(seqs, wt_seq, free_positions, bias_aas):
    """
    Detect native residues at free positions penalized by --bias_AAs.
    E.g., bias_AAs="E:-1.0" penalizes native E residues.
    """
    if not bias_aas:
        return []

    penalties = {}
    for token in bias_aas.split(","):
        token = token.strip()
        if ":" in token:
            aa, val = token.split(":")
            val = float(val)
            if val < 0:
                penalties[aa.upper()] = val

    if not penalties:
        return []

    victims = []
    for pos in free_positions:
        wt_aa = wt_seq[pos]
        if wt_aa.upper() not in penalties:
            continue

        design_aas = [s[pos] for s in seqs]
        counter = Counter(design_aas)
        n_wt = counter.get(wt_aa, 0)
        n_total = len(seqs)

        victims.append({
            "position_0idx": pos,
            "position_1idx": pos + 1,
            "wt_aa": wt_aa,
            "penalty": penalties[wt_aa.upper()],
            "n_designs_with_wt": n_wt,
            "n_total": n_total,
            "wt_fraction": n_wt / n_total,
            "distribution": dict(counter.most_common()),
            "forced_departure": n_wt == 0,
        })

    return victims


def compute_mi_matrix(seqs, positions):
    """
    Compute mutual information matrix with APC for given positions.
    Returns (mi_raw, mi_apc, position_list).
    """
    n_seqs = len(seqs)
    n_pos = len(positions)

    aa_to_idx = {aa: i for i, aa in enumerate("ACDEFGHIKLMNPQRSTVWY")}
    n_aa = 20

    freq_single = np.zeros((n_pos, n_aa))
    for seq in seqs:
        for pi, pos in enumerate(positions):
            aa = seq[pos]
            if aa in aa_to_idx:
                freq_single[pi, aa_to_idx[aa]] += 1
    freq_single /= n_seqs

    mi_raw = np.zeros((n_pos, n_pos))

    for i in range(n_pos):
        for j in range(i + 1, n_pos):
            freq_joint = np.zeros((n_aa, n_aa))
            for seq in seqs:
                aa_i = seq[positions[i]]
                aa_j = seq[positions[j]]
                if aa_i in aa_to_idx and aa_j in aa_to_idx:
                    freq_joint[aa_to_idx[aa_i], aa_to_idx[aa_j]] += 1
            freq_joint /= n_seqs

            mi = 0.0
            for a in range(n_aa):
                for b in range(n_aa):
                    if freq_joint[a, b] > 0 and freq_single[i, a] > 0 and freq_single[j, b] > 0:
                        mi += freq_joint[a, b] * math.log2(
                            freq_joint[a, b] / (freq_single[i, a] * freq_single[j, b])
                        )
            mi_raw[i, j] = mi
            mi_raw[j, i] = mi

    # APC correction
    mi_apc = np.zeros_like(mi_raw)
    row_means = mi_raw.mean(axis=1)
    global_mean = mi_raw.mean()

    if global_mean > 0:
        for i in range(n_pos):
            for j in range(i + 1, n_pos):
                apc = row_means[i] * row_means[j] / global_mean
                mi_apc[i, j] = mi_raw[i, j] - apc
                mi_apc[j, i] = mi_apc[i, j]

    return mi_raw, mi_apc, positions


def detect_artificial_couplings(seqs, wt_seq, _free_positions, victims, mi_apc,
                                 positions, mi_threshold=0.4):
    """
    Find high-MI pairs involving victim positions — these are likely
    artificial compensatory couplings created by the constraint.
    """
    if not victims:
        return []

    victim_pos_set = set(v["position_0idx"] for v in victims)
    pos_to_idx = {p: i for i, p in enumerate(positions)}

    couplings = []
    for v in victims:
        vpos = v["position_0idx"]
        if vpos not in pos_to_idx:
            continue
        vi = pos_to_idx[vpos]

        for j, pos_j in enumerate(positions):
            if j == vi:
                continue
            mi_val = mi_apc[vi, j]
            if mi_val >= mi_threshold:
                design_aas_j = [s[pos_j] for s in seqs]
                counter_j = Counter(design_aas_j)
                wt_aa_j = wt_seq[pos_j] if pos_j < len(wt_seq) else "?"

                couplings.append({
                    "victim_pos_1idx": vpos + 1,
                    "victim_wt_aa": v["wt_aa"],
                    "partner_pos_1idx": pos_j + 1,
                    "partner_wt_aa": wt_aa_j,
                    "mi_apc": round(float(mi_val), 4),
                    "partner_distribution": dict(counter_j.most_common(5)),
                    "partner_also_victim": pos_j in victim_pos_set,
                    "covariation_pairs": _extract_covariation_pairs(
                        seqs, vpos, pos_j
                    ),
                })

    couplings.sort(key=lambda x: -x["mi_apc"])
    return couplings


def _extract_covariation_pairs(seqs, pos_i, pos_j):
    """Show which AA pairs always co-occur."""
    pair_counts = Counter()
    for s in seqs:
        pair_counts[(s[pos_i], s[pos_j])] += 1

    result = []
    for (aa_i, aa_j), count in pair_counts.most_common(6):
        result.append({
            "pair": f"{aa_i}{pos_i+1}-{aa_j}{pos_j+1}",
            "count": count,
            "fraction": round(count / len(seqs), 3),
        })
    return result


def recommend_fixes(victims, couplings, mi_threshold=0.4):
    """
    Generate recommendations for which positions to add to fixed set.

    Severity logic:
      CRITICAL — high MI coupling exists (≥0.8), regardless of forced/reduced
      HIGH     — moderate MI coupling (≥threshold) OR forced departure
      MODERATE — forced departure but no coupling
      LOW      — reduced but not forced, no strong coupling
      UNAFFECTED — constraint had no actual effect (WT fraction = 100%)
    """
    recommendations = []

    for v in victims:
        pos = v["position_1idx"]
        wt = v["wt_aa"]

        if v["wt_fraction"] >= 1.0:
            continue

        involved_couplings = [
            c for c in couplings if c["victim_pos_1idx"] == pos
        ]
        max_mi = max((c["mi_apc"] for c in involved_couplings), default=0)

        if max_mi >= 0.8:
            severity = "CRITICAL"
        elif max_mi >= mi_threshold or (v["forced_departure"] and max_mi >= 0.3):
            severity = "HIGH"
        elif v["forced_departure"]:
            severity = "MODERATE"
        elif v["wt_fraction"] < 0.5:
            severity = "LOW"
        else:
            continue

        rec = {
            "action": f"FIX position {pos} to {wt}",
            "severity": severity,
            "reason": f"Native {wt}{pos} is {'completely excluded' if v['forced_departure'] else 'penalized'} "
                      f"by constraint. ",
            "n_artificial_couplings": len(involved_couplings),
            "max_mi_apc": round(max_mi, 4),
            "partner_positions": [c["partner_pos_1idx"] for c in involved_couplings],
        }

        if involved_couplings:
            rec["reason"] += (
                f"Creates {len(involved_couplings)} artificial coupling(s) "
                f"(max MI_APC={max_mi:.3f}) with positions "
                f"{', '.join(str(c['partner_pos_1idx']) for c in involved_couplings[:5])}."
            )
        else:
            rec["reason"] += "No strong artificial couplings detected, but forced departure from native is still a risk."

        recommendations.append(rec)

    # Also flag partner positions that could optionally be fixed
    victim_positions = set(v["position_1idx"] for v in victims)
    partner_fixes = {}
    for c in couplings:
        ppos = c["partner_pos_1idx"]
        if ppos not in victim_positions and ppos not in partner_fixes:
            partner_fixes[ppos] = {
                "action": f"OPTIONALLY FIX position {ppos} to {c['partner_wt_aa']}",
                "severity": "OPTIONAL",
                "reason": f"Coupled to victim position {c['victim_pos_1idx']} "
                          f"(MI_APC={c['mi_apc']:.3f}). Fixing the victim may "
                          f"resolve this coupling automatically.",
                "n_artificial_couplings": 0,
                "max_mi_apc": c["mi_apc"],
                "partner_positions": [c["victim_pos_1idx"]],
            }

    recommendations.extend(partner_fixes.values())
    recommendations.sort(key=lambda x: ("CRITICAL", "HIGH", "MODERATE", "LOW", "OPTIONAL").index(x["severity"]))

    return recommendations


def parse_msa(msa_path, wt_seq):
    """
    Parse a FASTA MSA and build a mapping from WT sequence positions (0-indexed)
    to alignment columns. Returns (aligned_seqs, wt_to_col mapping).

    The WT sequence must be the first entry in the MSA, OR the entry whose
    ungapped sequence best matches the provided wt_seq.
    """
    seqs = []
    headers = []
    current_header = None
    current_seq = []
    with open(msa_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header is not None:
                    headers.append(current_header)
                    seqs.append("".join(current_seq))
                current_header = line[1:]
                current_seq = []
            elif current_header is not None:
                current_seq.append(line)
        if current_header is not None:
            headers.append(current_header)
            seqs.append("".join(current_seq))

    if not seqs:
        raise ValueError(f"No sequences found in {msa_path}")

    # Try exact WT match first
    for i, s in enumerate(seqs):
        ungapped = s.replace("-", "").replace(".", "")
        if ungapped == wt_seq:
            wt_aligned = s
            wt_to_col = {}
            wt_pos = 0
            for col_idx, aa in enumerate(wt_aligned):
                if aa not in ("-", "."):
                    wt_to_col[wt_pos] = col_idx
                    wt_pos += 1
            print(f"  MSA: {len(seqs)} sequences, {len(wt_aligned)} columns, "
                  f"WT exact match: {headers[i][:40]}",
                  file=sys.stderr)
            return seqs, wt_to_col, headers

    # No exact match — find best (entry, offset) pair
    # Phase 1: coarse offset search (step=5) on ALL entries
    best_score = 0
    best_offset = 0
    best_match_idx = 0
    for i, s in enumerate(seqs):
        ungapped = s.replace("-", "").replace(".", "")
        ug_len = len(ungapped)
        max_off = max(1, len(wt_seq) - ug_len + 1)
        for off in range(0, max_off, 5):
            score = sum(1 for a, b in zip(wt_seq[off:off + ug_len], ungapped) if a == b)
            if score > best_score:
                best_score = score
                best_offset = off
                best_match_idx = i

    # Phase 2: fine-tune offset ±5 around best on the winning entry
    s = seqs[best_match_idx]
    ungapped = s.replace("-", "").replace(".", "")
    ug_len = len(ungapped)
    for off in range(max(0, best_offset - 5), min(len(wt_seq) - ug_len + 1, best_offset + 6)):
        score = sum(1 for a, b in zip(wt_seq[off:off + ug_len], ungapped) if a == b)
        if score > best_score:
            best_score = score
            best_offset = off

    best_match_score = best_score

    wt_aligned = seqs[best_match_idx]
    match_ungapped = wt_aligned.replace("-", "").replace(".", "")

    print(f"  MSA: {len(seqs)} sequences, {len(wt_aligned)} columns, "
          f"WT match: {headers[best_match_idx][:40]} "
          f"({best_match_score}/{len(match_ungapped)} identity, wt_offset={best_offset})",
          file=sys.stderr)

    wt_to_col = {}
    match_pos = 0
    for col_idx, aa in enumerate(wt_aligned):
        if aa not in ("-", "."):
            wt_pos = best_offset + match_pos
            if wt_pos < len(wt_seq):
                wt_to_col[wt_pos] = col_idx
            match_pos += 1

    return seqs, wt_to_col, headers


def evolutionary_validation(couplings, victims, msa_seqs, wt_to_col):
    """
    For each coupling and each victim, check evolutionary support
    by looking at amino acid frequencies in the natural MSA.

    For couplings: check if the MPNN-preferred AA pairs occur together
    in natural sequences.

    For victims: check if the MPNN-preferred substitution occurs in nature.

    Returns enriched couplings and victims with evolutionary annotations.
    """
    gap_chars = set("-.")

    def col_distribution(col_idx):
        counts = Counter()
        for s in msa_seqs:
            if col_idx < len(s) and s[col_idx] not in gap_chars:
                counts[s[col_idx]] += 1
        return counts

    def pair_distribution(col_i, col_j):
        counts = Counter()
        for s in msa_seqs:
            if (col_i < len(s) and col_j < len(s)
                    and s[col_i] not in gap_chars and s[col_j] not in gap_chars):
                counts[(s[col_i], s[col_j])] += 1
        return counts

    def classify_support(freq, total):
        if total == 0:
            return "NO_DATA"
        pct = freq / total * 100
        if pct < 0.1:
            return "NONE"
        if pct < 1.0:
            return "TRACE"
        if pct < 5.0:
            return "WEAK"
        if pct < 20.0:
            return "MODERATE"
        return "STRONG"

    # Validate each coupling
    for c in couplings:
        pos_i_0 = c["victim_pos_1idx"] - 1
        pos_j_0 = c["partner_pos_1idx"] - 1

        if pos_i_0 not in wt_to_col or pos_j_0 not in wt_to_col:
            c["evo_validation"] = {"status": "UNMAPPED", "detail": "Position not in MSA"}
            continue

        col_i = wt_to_col[pos_i_0]
        col_j = wt_to_col[pos_j_0]

        pair_dist = pair_distribution(col_i, col_j)
        total_pairs = sum(pair_dist.values())

        pair_annotations = []
        for cpair in c["covariation_pairs"]:
            pair_str = cpair["pair"]
            aa_i = pair_str[0]
            aa_j = pair_str.split("-")[1][0]
            natural_count = pair_dist.get((aa_i, aa_j), 0)

            pair_annotations.append({
                "mpnn_pair": f"{aa_i}-{aa_j}",
                "mpnn_fraction": cpair["fraction"],
                "natural_count": natural_count,
                "natural_total": total_pairs,
                "natural_pct": round(natural_count / total_pairs * 100, 2) if total_pairs > 0 else 0,
                "support": classify_support(natural_count, total_pairs),
            })

        dominant_pair = pair_annotations[0] if pair_annotations else None
        if dominant_pair:
            c["evo_validation"] = {
                "status": "CONFIRMED_ARTIFICIAL" if dominant_pair["support"] in ("NONE", "TRACE") else
                          "LIKELY_ARTIFICIAL" if dominant_pair["support"] == "WEAK" else
                          "POSSIBLY_GENUINE",
                "dominant_pair_natural_pct": dominant_pair["natural_pct"],
                "detail": (f"Dominant MPNN pair {dominant_pair['mpnn_pair']} "
                           f"({dominant_pair['mpnn_fraction']:.0%} of designs) found in "
                           f"{dominant_pair['natural_count']}/{dominant_pair['natural_total']} "
                           f"natural sequences ({dominant_pair['natural_pct']:.1f}%)"),
                "pair_annotations": pair_annotations,
            }
        else:
            c["evo_validation"] = {"status": "NO_DATA"}

    # Validate each victim's substitution pattern
    for v in victims:
        pos_0 = v["position_0idx"]
        if pos_0 not in wt_to_col:
            v["evo_validation"] = {"status": "UNMAPPED"}
            continue

        col = wt_to_col[pos_0]
        nat_dist = col_distribution(col)
        total = sum(nat_dist.values())

        top_mpnn = list(v["distribution"].keys())[0] if v["distribution"] else None
        if top_mpnn and top_mpnn != v["wt_aa"]:
            nat_count = nat_dist.get(top_mpnn, 0)
            v["evo_validation"] = {
                "wt_in_nature": nat_dist.get(v["wt_aa"], 0),
                "wt_natural_pct": round(nat_dist.get(v["wt_aa"], 0) / total * 100, 1) if total > 0 else 0,
                "mpnn_top_pick": top_mpnn,
                "mpnn_pick_in_nature": nat_count,
                "mpnn_pick_natural_pct": round(nat_count / total * 100, 2) if total > 0 else 0,
                "support_for_substitution": classify_support(nat_count, total),
                "natural_total": total,
                "natural_distribution": dict(nat_dist.most_common(5)),
            }
        else:
            v["evo_validation"] = {
                "wt_in_nature": nat_dist.get(v["wt_aa"], 0),
                "wt_natural_pct": round(nat_dist.get(v["wt_aa"], 0) / total * 100, 1) if total > 0 else 0,
                "natural_total": total,
                "natural_distribution": dict(nat_dist.most_common(5)),
            }


def compare_with_without_constraint(seqs_with, seqs_without, wt_seq, free_positions,
                                     omit_aas=None, bias_aas=None):
    """
    If both constrained and unconstrained checkpoints are available,
    compare to confirm which couplings are constraint-induced artifacts.
    """
    if omit_aas:
        omit_set = set(omit_aas.upper())
    else:
        omit_set = set()

    confirmed_artifacts = []

    for pos in free_positions:
        wt_aa = wt_seq[pos]
        if wt_aa.upper() not in omit_set:
            continue

        with_counter = Counter(s[pos] for s in seqs_with)
        without_counter = Counter(s[pos] for s in seqs_without)

        wt_frac_with = with_counter.get(wt_aa, 0) / len(seqs_with)
        wt_frac_without = without_counter.get(wt_aa, 0) / len(seqs_without)

        if wt_frac_with < wt_frac_without:
            confirmed_artifacts.append({
                "position_1idx": pos + 1,
                "wt_aa": wt_aa,
                "wt_fraction_constrained": round(wt_frac_with, 3),
                "wt_fraction_unconstrained": round(wt_frac_without, 3),
                "delta": round(wt_frac_without - wt_frac_with, 3),
                "distribution_constrained": dict(with_counter.most_common(5)),
                "distribution_unconstrained": dict(without_counter.most_common(5)),
            })

    confirmed_artifacts.sort(key=lambda x: -x["delta"])
    return confirmed_artifacts


def pre_design_scan(wt_seq, fixed_positions, omit_aas=None, bias_aas=None):
    """
    PRE-DESIGN MODE: Before running MPNN, predict which positions will
    create artificial correlations due to constraints. No ensemble needed.

    Returns dict with per-position warnings AND predicted correlation groups.
    """
    # fixed_positions are 1-indexed (from the scheme); convert to 0-indexed
    # to match the 0-indexed iteration below (free residues are reported as pos+1).
    fixed_set = set(p - 1 for p in fixed_positions) if fixed_positions else set()
    free = [i for i in range(len(wt_seq)) if i not in fixed_set]

    warnings = []
    groups = {}  # constraint_key -> [positions]

    if omit_aas:
        omit_set = set(omit_aas.upper())
        for pos in free:
            if wt_seq[pos].upper() in omit_set:
                aa = wt_seq[pos].upper()
                key = f"omit_{aa}"
                groups.setdefault(key, []).append(pos + 1)
                warnings.append({
                    "position_1idx": pos + 1,
                    "wt_aa": wt_seq[pos],
                    "constraint": f"omit_AAs={omit_aas}",
                    "impact": "FORCED_DEPARTURE",
                    "group": key,
                    "recommendation": (
                        f"Native {wt_seq[pos]}{pos+1} will be completely excluded. "
                        f"Add to fixed set OR accept non-native substitution."
                    ),
                })

    if bias_aas:
        for token in bias_aas.split(","):
            token = token.strip()
            if ":" in token:
                aa, val = token.split(":")
                val = float(val)
                if val < 0:
                    key = f"bias_{aa}_{val}"
                    for pos in free:
                        if wt_seq[pos].upper() == aa.upper():
                            groups.setdefault(key, []).append(pos + 1)
                            warnings.append({
                                "position_1idx": pos + 1,
                                "wt_aa": wt_seq[pos],
                                "constraint": f"bias_AAs {aa}:{val}",
                                "impact": "PENALIZED",
                                "group": key,
                                "recommendation": (
                                    f"Native {wt_seq[pos]}{pos+1} is penalized by {val}. "
                                    f"MPNN may substitute it. Consider fixing if essential."
                                ),
                            })

    # Predict correlation structure
    correlation_groups = []
    total_predicted_correlations = 0
    for key, positions in groups.items():
        n = len(positions)
        n_corr = n * (n - 1) // 2
        total_predicted_correlations += n_corr
        correlation_groups.append({
            "constraint": key,
            "n_positions": n,
            "positions": positions,
            "predicted_artificial_correlations": n_corr,
            "severity": "HIGH" if n >= 5 else ("MEDIUM" if n >= 3 else "LOW"),
        })

    return {
        "warnings": warnings,
        "correlation_groups": correlation_groups,
        "total_artifact_positions": len(warnings),
        "total_predicted_correlations": total_predicted_correlations,
    }


def run_analysis(checkpoint_path, wt_seq, fixed_positions=None,
                 omit_aas=None, bias_aas=None, comparison_checkpoint=None,
                 mi_threshold=0.4, msa_path=None):
    """Main analysis pipeline."""

    seqs = parse_checkpoint(checkpoint_path)
    seq_len = len(seqs[0])

    if len(wt_seq) != seq_len:
        print(f"WARNING: WT length ({len(wt_seq)}) != design length ({seq_len}). "
              f"Truncating/padding as needed.", file=sys.stderr)
        if len(wt_seq) > seq_len:
            wt_seq = wt_seq[:seq_len]
        else:
            wt_seq = wt_seq + "X" * (seq_len - len(wt_seq))

    fixed, free = identify_fixed_free(seqs, fixed_positions)

    report = {
        "checkpoint": str(checkpoint_path),
        "n_designs": len(seqs),
        "seq_length": seq_len,
        "n_fixed": len(fixed),
        "n_free": len(free),
        "constraints": {
            "omit_AAs": omit_aas,
            "bias_AAs": bias_aas,
        },
    }

    # Step 1: Detect omit victims
    omit_victims = detect_omit_victims(seqs, wt_seq, free, omit_aas)
    bias_victims = detect_bias_victims(seqs, wt_seq, free, bias_aas)
    all_victims = omit_victims + bias_victims

    report["omit_victims"] = omit_victims
    report["bias_victims"] = bias_victims
    report["n_omit_victims"] = len(omit_victims)
    report["n_bias_victims"] = len(bias_victims)

    # Step 2: Compute MI matrix on free positions
    print(f"Computing MI matrix for {len(free)} free positions across {len(seqs)} designs...",
          file=sys.stderr)
    _mi_raw, mi_apc, positions = compute_mi_matrix(seqs, free)

    # Step 3: Detect artificial couplings
    couplings = detect_artificial_couplings(
        seqs, wt_seq, free, all_victims, mi_apc, positions, mi_threshold
    )
    report["artificial_couplings"] = couplings
    report["n_artificial_couplings"] = len(couplings)

    # Step 4: Generate recommendations
    recommendations = recommend_fixes(all_victims, couplings, mi_threshold)
    report["recommendations"] = recommendations

    # Step 5: Evolutionary validation via MSA
    if msa_path:
        print(f"Running evolutionary validation from MSA...", file=sys.stderr)
        msa_seqs, wt_to_col, _ = parse_msa(msa_path, wt_seq)
        evolutionary_validation(couplings, all_victims, msa_seqs, wt_to_col)
        report["msa_path"] = str(msa_path)
        report["msa_n_sequences"] = len(msa_seqs)

        for r in recommendations:
            pos_1idx = int(r["action"].split("position ")[1].split(" ")[0])
            matching_couplings = [c for c in couplings
                                  if c.get("victim_pos_1idx") == pos_1idx
                                  and "evo_validation" in c]
            if matching_couplings:
                evo = matching_couplings[0]["evo_validation"]
                status = evo.get("status", "")
                if status == "CONFIRMED_ARTIFICIAL" and r["severity"] != "CRITICAL":
                    r["severity"] = "CRITICAL"
                    r["reason"] += f" Evolutionary validation: {status}."
                elif status in ("CONFIRMED_ARTIFICIAL", "LIKELY_ARTIFICIAL"):
                    r["reason"] += f" Evolutionary validation: {status}."

    # Step 6: Optional comparison with unconstrained checkpoint
    if comparison_checkpoint:
        seqs_without = parse_checkpoint(comparison_checkpoint)
        confirmed = compare_with_without_constraint(
            seqs, seqs_without, wt_seq, free, omit_aas, bias_aas
        )
        report["confirmed_artifacts"] = confirmed

    # Step 6: Global MI statistics for context
    flat_mi = mi_apc[np.triu_indices_from(mi_apc, k=1)]
    report["mi_stats"] = {
        "mean": round(float(flat_mi.mean()), 4),
        "std": round(float(flat_mi.std()), 4),
        "max": round(float(flat_mi.max()), 4),
        "n_pairs_above_threshold": int((flat_mi >= mi_threshold).sum()),
        "threshold": mi_threshold,
    }

    return report


def format_report(report):
    """Pretty-print the analysis report."""
    lines = []
    lines.append("=" * 72)
    lines.append("MPNN BLIND SPOT DETECTOR — ANALYSIS REPORT")
    lines.append("=" * 72)
    lines.append(f"Checkpoint: {report['checkpoint']}")
    lines.append(f"Designs: {report['n_designs']}  |  Length: {report['seq_length']}aa  |  "
                 f"Fixed: {report['n_fixed']}  |  Free: {report['n_free']}")
    lines.append(f"Constraints: omit_AAs={report['constraints']['omit_AAs']}  "
                 f"bias_AAs={report['constraints']['bias_AAs']}")
    lines.append("")

    # Victims (filter out positions where constraint had no effect)
    all_victims = [v for v in report["omit_victims"] + report["bias_victims"]
                   if v["wt_fraction"] < 1.0]
    unaffected = [v for v in report["omit_victims"] + report["bias_victims"]
                  if v["wt_fraction"] >= 1.0]
    if all_victims:
        lines.append(f"{'─'*72}")
        lines.append(f"NATIVE RESIDUES AFFECTED BY CONSTRAINTS ({len(all_victims)} found"
                     f"{f', {len(unaffected)} unaffected filtered out' if unaffected else ''})")
        lines.append(f"{'─'*72}")
        for v in sorted(all_victims, key=lambda x: x["wt_fraction"]):
            flag = "FORCED OUT" if v["forced_departure"] else "REDUCED"
            constraint = "omit_AAs" if v.get("omitted") else f"bias_AAs ({v.get('penalty', '?')})"
            top3 = ", ".join(f"{aa}:{n}" for aa, n in list(v["distribution"].items())[:4])
            lines.append(
                f"  [{flag}] {v['wt_aa']}{v['position_1idx']} "
                f"(WT fraction: {v['wt_fraction']:.0%}, constraint: {constraint})"
            )
            lines.append(f"           Distribution: {top3}")
    else:
        lines.append("No native residues affected by constraints.")

    # Artificial couplings
    if report["artificial_couplings"]:
        lines.append("")
        lines.append(f"{'─'*72}")
        lines.append(f"ARTIFICIAL CO-VARIATION PATTERNS ({report['n_artificial_couplings']} found)")
        lines.append(f"{'─'*72}")
        for c in report["artificial_couplings"][:10]:
            victim_tag = f"*VICTIM* {c['victim_wt_aa']}{c['victim_pos_1idx']}"
            partner_tag = f"{c['partner_wt_aa']}{c['partner_pos_1idx']}"
            if c["partner_also_victim"]:
                partner_tag += " *ALSO VICTIM*"

            evo_tag = ""
            if "evo_validation" in c:
                ev = c["evo_validation"]
                status = ev.get("status", "")
                if status in ("CONFIRMED_ARTIFICIAL", "LIKELY_ARTIFICIAL"):
                    evo_tag = f"  [{status}]"
                elif status == "POSSIBLY_GENUINE":
                    evo_tag = f"  [POSSIBLY GENUINE]"

            lines.append(f"  {victim_tag} <-> {partner_tag}  MI_APC={c['mi_apc']:.3f}{evo_tag}")
            for pair in c["covariation_pairs"][:4]:
                evo_note = ""
                if "evo_validation" in c and "pair_annotations" in c["evo_validation"]:
                    for pa in c["evo_validation"]["pair_annotations"]:
                        if pa["mpnn_pair"] == f"{pair['pair'][0]}-{pair['pair'].split('-')[1][0]}":
                            evo_note = f"  (natural: {pa['natural_pct']:.1f}% — {pa['support']})"
                            break
                lines.append(f"    {pair['pair']}: {pair['count']}/{report['n_designs']} "
                             f"({pair['fraction']:.0%}){evo_note}")

    # MI stats
    mi = report["mi_stats"]
    lines.append("")
    lines.append(f"MI stats: mean={mi['mean']:.4f}, std={mi['std']:.4f}, max={mi['max']:.4f}")
    lines.append(f"Pairs above threshold ({mi['threshold']}): {mi['n_pairs_above_threshold']}")

    # Confirmed artifacts (if comparison available)
    if "confirmed_artifacts" in report and report["confirmed_artifacts"]:
        lines.append("")
        lines.append(f"{'─'*72}")
        lines.append("CONFIRMED ARTIFACTS (constrained vs unconstrained comparison)")
        lines.append(f"{'─'*72}")
        for a in report["confirmed_artifacts"]:
            lines.append(
                f"  {a['wt_aa']}{a['position_1idx']}: "
                f"WT fraction {a['wt_fraction_unconstrained']:.0%} → "
                f"{a['wt_fraction_constrained']:.0%} (Δ={a['delta']:+.0%})"
            )

    # Recommendations
    if report["recommendations"]:
        lines.append("")
        lines.append(f"{'─'*72}")
        lines.append("RECOMMENDATIONS")
        lines.append(f"{'─'*72}")
        for r in report["recommendations"]:
            lines.append(f"  [{r['severity']}] {r['action']}")
            lines.append(f"    {r['reason']}")
            if r["partner_positions"]:
                lines.append(f"    Coupled to: {r['partner_positions']}")

    lines.append("")
    lines.append("=" * 72)
    return "\n".join(lines)


def format_pre_design_report(scan_result, wt_len, n_fixed, omit_aas, bias_aas):
    """Pretty-print pre-design scan results with correlation predictions."""
    lines = []
    lines.append("=" * 72)
    lines.append("MPNN BLIND SPOT DETECTOR — PRE-DESIGN SCAN")
    lines.append("=" * 72)
    lines.append(f"WT length: {wt_len}aa  |  Fixed: {n_fixed}  |  Free: {wt_len - n_fixed}")
    lines.append(f"Planned constraints: omit_AAs={omit_aas}  bias_AAs={bias_aas}")
    lines.append("")

    warnings = scan_result.get("warnings", [])
    groups = scan_result.get("correlation_groups", [])
    total_corr = scan_result.get("total_predicted_correlations", 0)

    forced = [w for w in warnings if w["impact"] == "FORCED_DEPARTURE"]
    penalized = [w for w in warnings if w["impact"] == "PENALIZED"]

    if forced:
        lines.append(f"{'─'*72}")
        lines.append(f"WILL BE FORCED OUT ({len(forced)} native residues)")
        lines.append(f"{'─'*72}")
        for w in forced:
            lines.append(f"  {w['wt_aa']}{w['position_1idx']} — {w['recommendation']}")

    if penalized:
        lines.append(f"{'─'*72}")
        lines.append(f"WILL BE PENALIZED ({len(penalized)} native residues)")
        lines.append(f"{'─'*72}")
        for w in penalized:
            lines.append(f"  {w['wt_aa']}{w['position_1idx']} — {w['recommendation']}")

    if groups:
        lines.append("")
        lines.append(f"{'─'*72}")
        lines.append(f"PREDICTED ARTIFICIAL CORRELATIONS ({total_corr} total)")
        lines.append(f"{'─'*72}")
        lines.append(f"  All positions sharing a constraint will co-vary in every")
        lines.append(f"  design, creating artificial LD/MI regardless of structure.")
        lines.append("")
        for g in sorted(groups, key=lambda x: -x["predicted_artificial_correlations"]):
            pos_str = ", ".join(str(p) for p in g["positions"][:8])
            if len(g["positions"]) > 8:
                pos_str += f", ... (+{len(g['positions'])-8} more)"
            lines.append(
                f"  [{g['severity']}] {g['constraint']}: "
                f"{g['n_positions']} positions -> "
                f"{g['predicted_artificial_correlations']} artificial pairs"
            )
            lines.append(f"         positions: {pos_str}")

    if not warnings:
        lines.append("No native residues at risk. Constraints are safe.")

    lines.append("")
    lines.append("=" * 72)
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Detect MPNN blind spots: native residues forcibly mutated by constraints"
    )
    parser.add_argument("checkpoint", nargs="?", help="MPNN checkpoint JSON file (omit for --pre-design mode)")
    parser.add_argument("--wt", help="Wild-type protein sequence (single string)")
    parser.add_argument("--wt_file", help="FASTA file containing WT sequence (first entry used)")
    parser.add_argument("--omit_AAs", help="Amino acids omitted during design (e.g., 'C')")
    parser.add_argument("--bias_AAs", help="Bias string used during design (e.g., 'E:-1.0')")
    parser.add_argument("--fixed", help="Comma-separated 0-indexed fixed positions")
    parser.add_argument("--scheme", help="Fixation scheme JSON file")
    parser.add_argument("--path", help="Path name within fixation scheme")
    parser.add_argument("--comparison", help="Unconstrained checkpoint for artifact confirmation")
    parser.add_argument("--mi_threshold", type=float, default=0.4,
                        help="MI_APC threshold for coupling detection (default: 0.4)")
    parser.add_argument("--output", help="Output JSON file path")
    parser.add_argument("--quiet", action="store_true", help="Suppress text report")
    parser.add_argument("--msa", help="FASTA MSA file for evolutionary validation of couplings")
    parser.add_argument("--pre_design", action="store_true",
                        help="Pre-design scan mode: check constraints BEFORE running MPNN")

    args = parser.parse_args()

    # Load WT sequence
    if args.wt_file:
        with open(args.wt_file) as f:
            wt_seq = ""
            for line in f:
                if not line.startswith(">"):
                    wt_seq += line.strip()
        if not wt_seq:
            print("ERROR: Could not parse WT sequence from FASTA", file=sys.stderr)
            sys.exit(1)
    else:
        wt_seq = args.wt

    # Load fixed positions
    fixed_positions = None
    if args.scheme and args.path:
        with open(args.scheme) as f:
            scheme = json.load(f)
        if args.path in scheme:
            fixed_positions = scheme[args.path]["fixed_positions"]
        else:
            print(f"ERROR: Path '{args.path}' not found in scheme. "
                  f"Available: {list(scheme.keys())}", file=sys.stderr)
            sys.exit(1)
    elif args.fixed:
        fixed_positions = [int(x) for x in args.fixed.split(",")]

    # Load constraints from scheme if available
    omit_aas = args.omit_AAs
    bias_aas = args.bias_AAs
    if args.scheme and args.path:
        with open(args.scheme) as f:
            scheme = json.load(f)
        path_data = scheme.get(args.path, {})
        omit_aas = omit_aas or path_data.get("omit_amino_acids")
        bias_aas = bias_aas or path_data.get("bias_amino_acids")

    if not wt_seq:
        print("ERROR: --wt or --wt_file required", file=sys.stderr)
        sys.exit(1)

    # Pre-design mode: no checkpoint needed
    if args.pre_design:
        scan_result = pre_design_scan(
            wt_seq,
            fixed_positions or [],
            omit_aas=omit_aas,
            bias_aas=bias_aas,
        )
        if not args.quiet:
            n_fixed = len(fixed_positions) if fixed_positions else 0
            print(format_pre_design_report(scan_result, len(wt_seq), n_fixed, omit_aas, bias_aas))
        if args.output:
            with open(args.output, "w") as f:
                json.dump(scan_result, f, indent=2)
            print(f"\nJSON report saved to: {args.output}", file=sys.stderr)
        return scan_result

    if not args.checkpoint:
        print("ERROR: checkpoint required (use --pre_design for scan without checkpoint)",
              file=sys.stderr)
        sys.exit(1)

    report = run_analysis(
        args.checkpoint,
        wt_seq,
        fixed_positions=fixed_positions,
        omit_aas=omit_aas,
        bias_aas=bias_aas,
        comparison_checkpoint=args.comparison,
        mi_threshold=args.mi_threshold,
        msa_path=args.msa,
    )

    if not args.quiet:
        print(format_report(report))

    if args.output:
        with open(args.output, "w") as f:
            json.dump(report, f, indent=2)
        print(f"\nJSON report saved to: {args.output}", file=sys.stderr)

    return report


if __name__ == "__main__":
    main()
