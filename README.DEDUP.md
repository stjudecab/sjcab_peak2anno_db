# Deduplicate Selection methods:

- `longcol5`: no selector is needed; select the isoform with the largest
  numeric value in BED column 5. This is the default for `dedup-bed`.
- `long`: no selector is needed; select the isoform with the largest `end - start`.
- `peak`: selector is a peak BED file with peak score in column 5; the isoform
  whose TSS +/- promoter window has the highest peak score is selected. Text
  selectors are also accepted, for example `chr1:100-200 10`.
- `isoID`: selector is a transcript ID list, one ID per line.
- `isoexp`: selector is a two-column table: transcript ID, then expression.
- `perover`: selector is a BED file, often user-filtered ChromHMM active states;
  the isoform promoter with the highest percent overlap is selected. Text
  selectors are also accepted, for example `chr1:100-200` or `chr1_100_200`.

# Defaults and naming:

- Promoter half-window for `peak` and `perover` is `2kb`; change with
  `--promoter-bp`/`-p`.
- Matching is inclusive by default: versioned and unversioned transcript IDs can
  match each other, and any BED overlap counts.
- `--exclusive` requires exact transcript ID matches for `isoID`/`isoexp` and
  BED features fully contained inside the promoter for `peak`/`perover`.
- Text `peak`/`perover` selectors may have a header or no header. Region strings
  accept the common `sjcab_peak2anno` delimiters, including `:`, `-`, `_`, `/`,
  `;`, and `,`.
- Isoforms are grouped by gene symbol by default; use `--gene-key ensid` to
  group by Ensembl/GENCODE gene ID.
- Species/version lookup writes prefixes such as
  `hg38.v31.deduppeak` or `hg38.v31.filterperover`.
- Explicit BED input such as `my.bed` writes prefixes such as
  `my.deduppeak` or `my.filterperover`.

# Python API:

```python
import sjcab_peak2anno_db as db

db.dedup_bed(
    "peak",
    "h3k4me3_peaks.bed",
    output_dir="annotations",
    species="hg38",
    version="v31",
)

db.filter_bed(
    "perover",
    "active_chromhmm.bed",
    output_dir="annotations",
    gene_bed="annotations/bed/hg38/v31/all.gene.bed",
    promoter_bp="2kb",
)
```

# Deduplicate Selector Examples

These examples use an all-isoform BED. BED column 4 is the gene name used to
group transcripts. `GENE1` intentionally has two non-overlapping transcript
regions so selector results are visible.

Input gene BED:

```text
chr1    100    200    GENE1    100    +    ENSG1    ENST1
chr1    500    600    GENE1    100    +    ENSG1    ENST2
chr1    800    900    GENE2    100    -    ENSG2    ENST3
```

## `longcol5`

The default `dedup-bed` selector. It keeps the row with the largest numeric
value in BED column 5 for each gene.

```text
chr1    100    200    GENE1    100    +    ENSG1    ENST1
chr1    500    600    GENE1    250    +    ENSG1    ENST2
chr1    800    900    GENE2    100    -    ENSG2    ENST3
```

```bash
# Keep the largest value in BED column 5 for each gene.
sjcab-peak2anno-db dedup-bed hg38 v31 -m longcol5 -b genes.bed -o .
```

Output:

```text
chr1    500    600    GENE1    250    +    ENSG1    ENST2
chr1    800    900    GENE2    100    -    ENSG2    ENST3
```

## `long`

Keeps the row with the largest interval length, `end - start`, ignoring
column 5. These intervals tie, so input order resolves the tie.

```text
chr1    100    200    GENE1    20     +    ENSG1    ENST1
chr1    500    600    GENE1    999    +    ENSG1    ENST2
chr1    800    900    GENE2    100    -    ENSG2    ENST3
```

```bash
# Keep the longest genomic interval for each gene.
sjcab-peak2anno-db dedup-bed hg38 v31 -m long -b genes.bed -o .
```

Output:

```text
chr1    100    200    GENE1    20     +    ENSG1    ENST1
chr1    800    900    GENE2    100    -    ENSG2    ENST3
```

## `peak`

The peak score is BED column 5. Both `GENE1` regions are non-overlapping, but
the second region overlaps the higher-score peak, so `ENST2` is kept.

Selector BED before selection:

```text
chr1    105    115    peak1    10
chr1    505    515    peak2    50
```

The same selector may be whitespace-delimited text, with or without a header:

```text
region score
chr1:105-115 10
chr1_505_515 50
```

```bash
# Select the transcript whose promoter overlaps the highest-score peak.
sjcab-peak2anno-db dedup-bed hg38 v31 -m peak -p 20 -b genes.bed -i peaks.txt -o .
```

Output:

```text
chr1    500    600    GENE1    100    +    ENSG1    ENST2
chr1    800    900    GENE2    100    -    ENSG2    ENST3
```

## `perover`

The selector score is percent overlap. The second `GENE1` region has the
larger overlap, so `ENST2` is kept.

Selector BED before selection:

```text
chr1    100    110    active1
chr1    500    590    active2
```

Headered or headerless text is also accepted:

```text
region
chr1:100-110
chr1/500/590
```

```bash
# Select the transcript with the greatest selector overlap.
sjcab-peak2anno-db filter-bed -m perover -p 20 -b genes.bed -i active.txt -o .
```

Output:

```text
chr1    500    600    GENE1    100    +    ENSG1    ENST2
```

## `isoID`

The selector is a transcript ID list. `dedup-bed` falls back to the longest
row for genes without a match.

```text
ENST2
ENST3
```

```bash
# Keep the transcript IDs listed in ids.txt.
sjcab-peak2anno-db dedup-bed hg38 v31 -m isoID -b genes.bed -i ids.txt -o .
```

Output:

```text
chr1    500    600    GENE1    100    +    ENSG1    ENST2
chr1    800    900    GENE2    100    -    ENSG2    ENST3
```

## `isoexp`

The selector is a two-column transcript ID and expression table. The first row
is a header and is ignored.

```text
transcript_id    expression
ENST1            2.0
ENST2            8.5
ENST3            1.0
```

```bash
# Keep the highest-expression transcript for each gene.
sjcab-peak2anno-db dedup-bed hg38 v31 -m isoexp -b genes.bed -i expression.tsv -o .
```

Output:

```text
chr1    500    600    GENE1    100    +    ENSG1    ENST2
chr1    800    900    GENE2    100    -    ENSG2    ENST3
```

Region text accepts the common `sjcab_peak2anno` delimiters, including `:`,
`-`, `_`, `/`, `;`, and `,`.
