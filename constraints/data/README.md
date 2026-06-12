# constraint input data

## `cas9_msa.fasta.gz`
Cas9 family multiple sequence alignment used for both the conservation and the
coevolution/DCA stages.

- **3,782 sequences**, alignment length 1,897
- **Row 0 = `SpCas9_WT_query`** — the reference. Its ungapped length (1,341)
  matches the 8SRS chain-A modelled residues (author numbering 3–1366), which is
  what makes the gap-safe column→PDB mapping exact (0 non-SpRY mismatches).
- gzip'd (388 KB vs 7 MB) — `msa_io.parse_msa` reads `.gz` transparently.

Both signals are derived from this single alignment, so conservation (1st-order)
and DCA (2nd-order) are strictly apples-to-apples.
