"""
Cas9 design-time constraint generation.

Permanent pipeline stage that turns a family MSA + a target structure into the
fixed-position lists fed to LigandMPNN. Two complementary signals:

  conservation  (1st-order, single-column)  -> conservation.py
  coevolution   (2nd-order, column pairs)    -> coevolution.py   [DCA / EVcouplings]

Both are mapped to PDB author numbering by mapping.py and merged into tiers by
pipeline.py. This module formalizes the previously ad-hoc /tmp DCA scripts
(Task #41).
"""
from .conservation import compute_conservation
from .coevolution import compute_coevolution
from .mapping import build_column_to_pdb

__all__ = ["compute_conservation", "compute_coevolution", "build_column_to_pdb"]
