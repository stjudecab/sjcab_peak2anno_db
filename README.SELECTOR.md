# Selector Examples

These examples use an all-isoform gene BED. Each gene may have multiple
transcripts; the selector chooses one transcript per gene.

## Length Selectors

`longcol5` is the default for `dedup-bed`. It reads the numeric value in BED
column 5, which is normally the transcript's exon length:

```text
chr1    100    300    TX1    180    +    GENE1
chr1    100    280    TX2    160    +    GENE1
```

```bash
sjcab-peak2anno-db dedup-bed hg38 v31 -b genes.bed
sjcab-peak2anno-db dedup-bed hg38 v31 -m longcol5 -b genes.bed
```

`long` ignores column 5 and selects the largest interval length, `end - start`:

```text
chr1    100    300    TX1    20     +    GENE1
chr1    100    280    TX2    999    +    GENE1
```

```bash
sjcab-peak2anno-db dedup-bed hg38 v31 -m long -b genes.bed
```

## Peak Selectors

BED input uses score column 5:

```text
chr1    120    130    peak1    10
chr1    240    260    peak2    25
```

Text input may have a header or no header. The first field is a region and the
next numeric field is its score:

```text
region score
chr1:120-130 10
chr1_240_260 25
```

```bash
sjcab-peak2anno-db dedup-bed hg38 v31 -m peak -b genes.bed -i peaks.txt
```

## Percent-Overlap Selectors

BED input:

```text
chr1    100    180    active1
chr1    240    270    active2
```

Headered or headerless text input:

```text
region
chr1:100-180
chr1/240/270
```

```bash
sjcab-peak2anno-db filter-bed -b genes.bed -m perover -i active.txt
```

Region text accepts the same common delimiters as `sjcab_peak2anno`, including
`:`, `-`, `_`, `/`, `;`, and `,`.
