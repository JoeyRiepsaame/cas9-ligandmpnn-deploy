"""
Gap-safe mapping between MSA columns and PDB author residue numbers.

This is the stage the original /tmp DCA script got WRONG: it labelled positions
1,2,3,... by ungapped MSA index, but 8SRS chain A starts at author residue 3 and
contains internal gaps. The fix: the i-th ungapped residue of the WT query maps
to the i-th *present* author residue number in the structure, which naturally
jumps across missing loops.

Canonical key throughout the module: the ORIGINAL MSA column index. A column maps
to a PDB residue iff the WT query (row 0) is non-gap there. Columns where WT is a
gap have no structural counterpart and are intentionally absent from the map.
"""
from __future__ import annotations

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E",
    "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}


def parse_structure_chain(path, chain="A"):
    """{author_resnum: one_letter_aa} for protein residues of `chain`. .pdb or .cif."""
    return _parse_cif(path, chain) if str(path).endswith(".cif") else _parse_pdb(path, chain)


def _parse_pdb(path, chain):
    res = {}
    with open(path) as f:
        for line in f:
            if line[:4] == "ATOM" and line[21] == chain:
                num = int(line[22:26])
                name = line[17:20].strip()
                if num not in res and name in THREE_TO_ONE:
                    res[num] = THREE_TO_ONE[name]
    return res


def _parse_cif(path, chain):
    res = {}
    in_atom, cols, idx = False, [], {}
    with open(path) as f:
        for line in f:
            if line.startswith("_atom_site."):
                cols.append(line.strip().split(".")[1]); in_atom = True; continue
            if in_atom and line.startswith("ATOM"):
                if not idx:
                    idx = {c: i for i, c in enumerate(cols)}
                p = line.split()
                if p[idx["auth_asym_id"]] != chain:
                    continue
                name = p[idx["auth_comp_id"]]
                if name not in THREE_TO_ONE:
                    continue
                num = int(p[idx["auth_seq_id"]])
                if num not in res:
                    res[num] = THREE_TO_ONE[name]
            elif in_atom and line.startswith("#"):
                in_atom = False
    return res


def build_column_to_pdb(wt_aligned_seq, structure_residues,
                        spry_positions=None, strict=True):
    """
    Map ORIGINAL MSA column index -> PDB author residue number (gap-safe).

    wt_aligned_seq    : aligned (gapped) WT query sequence (row 0 of the MSA)
    structure_residues: {author_resnum: aa} from parse_structure_chain
    spry_positions    : author resnums where WT(SpCas9) != structure(SpRY) is expected.

    Returns (col_to_pdb: dict[int col -> int author_resnum], report: dict).
    """
    spry = set(spry_positions or [])
    pdb_sorted = sorted(structure_residues)

    col_to_pdb, mismatches = {}, []
    ungapped_idx = 0
    for col, aa in enumerate(wt_aligned_seq):
        if aa == "-":
            continue
        if ungapped_idx >= len(pdb_sorted):
            raise ValueError("WT query has more residues than the structure chain")
        pdb_num = pdb_sorted[ungapped_idx]
        col_to_pdb[col] = pdb_num
        struct_aa = structure_residues[pdb_num]
        if struct_aa != aa and pdb_num not in spry:
            mismatches.append((ungapped_idx, pdb_num, struct_aa, aa))
        ungapped_idx += 1

    report = {
        "wt_ungapped_len": ungapped_idx,
        "structure_residues": len(pdb_sorted),
        "structure_range": (pdb_sorted[0], pdb_sorted[-1]) if pdb_sorted else None,
        "columns_mapped": len(col_to_pdb),
        "unexpected_mismatches": mismatches,
    }
    if strict and mismatches:
        raise ValueError(
            f"{len(mismatches)} non-SpRY WT/structure mismatches — MSA and "
            f"structure are misaligned. First few: {mismatches[:5]}"
        )
    return col_to_pdb, report
