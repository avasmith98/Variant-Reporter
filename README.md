# Variant Reporter

Run `python3 Variant_Reporter.py` with the dependencies in `requirements.txt` installed.

## Transcript and variant verification

- Input includes a versioned RefSeq transcript, gene, cDNA change and protein change, for example `NM_005228.3(EGFR):c.2648T>C (p.Leu883Ser)`.
- NCBI supplies the current version of that same transcript accession.
- ClinVar must contain the exact original cDNA expression and one expression for the current transcript under the same simple allele. Updated cDNA and protein annotations come from that record; coordinates are never copied onto a new transcript version.
- The original protein change must match its paired ClinVar annotation. When historical protein annotations are absent, coding substitutions can be checked against the original version's RefSeq codon. Other unverifiable changes fail.
- Ensembl must provide exactly one explicit MANE mapping for the current RefSeq accession **and version** on GRCh38. Non-MANE or missing mappings fail; the first search result is never used as a substitute.
- Genomic location and classifications come from the matched ClinVar record. Reports show the original input alongside the updated annotations.
- Missing required data, ambiguous mappings and service failures prevent PDF generation. Protein-domain requests retry twice, then fail. Missing domain descriptions also fail rather than becoming blank text.

This uses recorded equivalence from ClinVar, not a general-purpose transcript alignment/remapping algorithm. A valid variant that lacks the required source annotations will be rejected rather than guessed.

Run regression tests with `python3 -m unittest discover -s tests -v`.
The EGFR XML fixture contains selected fields retrieved from NCBI ClinVar variation 1050891 via EFetch on 2026-09-29. Tests that modify it represent synthetic cases.
