# README

Filter sequences for the human v2 RiboFlow filter pass. RiboFlow aligns reads
to these sequences before transcriptome alignment and discards anything that
matches, so unwanted RNA species and repeat-derived reads never reach the
transcriptome.

## Sequences added (relative to v1)

v2 extends the v1 filter with three small-RNA pseudogenes (U4, U6, and 7SK
families) and two RepeatMasker consensus elements (a LINE and a SINE), so
reads originating from these elements are stripped before transcriptome
alignment.

| Sequence ID | Locus (hg38) | Strand | Source / type |
| --- | --- | --- | --- |
| `hg38_RNU4-46Px_ENST00000410818.1` | chr16:19,498,610–19,498,747 | + | RNU4-46P — U4 snRNA pseudogene (Ensembl) |
| `hg38_RNU6-1330P_ENST00000362775.1` | chr5:74,779,309–74,779,413 | + | RNU6-1330P — U6 snRNA pseudogene (Ensembl) |
| `hg38_RN7SKP176_ENST00000363673.1` | chr16:81,961,926–81,962,243 | + | RN7SKP176 — 7SK RNA pseudogene (Ensembl) |
| `hg38_rmsk_L1MB5` | chr2:175,176,603–175,177,048 | − | RepeatMasker L1MB5 (LINE-1 subfamily) |
| `hg38_rmsk_AluY` | chr17:75,305,393–75,305,688 | − | RepeatMasker AluY (SINE / Alu) |