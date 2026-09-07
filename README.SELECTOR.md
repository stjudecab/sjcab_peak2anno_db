# Selector Examples

All examples use an all-isoform BED. BED column 4 is the gene name used to
group transcripts:

```text
chr1    100    300    GENE1    180    +    ENSG1    ENST1
chr1    100    280    GENE1    160    +    ENSG1    ENST2
chr1    500    620    GENE2    100    -    ENSG2    ENST3
```

## `longcol5`

The default selector. It selects the largest numeric value in BED column 5.

```bash
sjcab-peak2anno-db dedup-bed hg38 v31 -b genes.bed -o .
sjcab-peak2anno-db dedup-bed hg38 v31 -m longcol5 -b genes.bed -o .
```

Result:

```text
chr1    100    300    GENE1    180    +    ENSG1    ENST1
chr1    500    620    GENE2    100    -    ENSG2    ENST3
```

## `long`

Selects the largest interval length, `end - start`, ignoring column 5.

Input:

```text
chr1    100    300    GENE1    20     +    ENSG1    ENST1
chr1    100    280    GENE1    999    +    ENSG1    ENST2
```

```bash
sjcab-peak2anno-db dedup-bed hg38 v31 -m long -b genes.bed -o .
```

Result:

```text
chr1    100    300    GENE1    20     +    ENSG1    ENST1
```

## `peak`

BED selector input uses score column 5:

```text
chr1    120    130    peak1    10
chr1    240    260    peak2    25
```

Text selectors may have a header or no header. The first field is a region and
the next numeric field is its score:

```text
region score
chr1:120-130 10
chr1_240_260 25
```

```bash
sjcab-peak2anno-db dedup-bed hg38 v31 -m peak -b genes.bed -i peaks.txt -o .
```

Result:

```text
chr1    100    300    GENE1    180    +    ENSG1    ENST1
chr1    500    620    GENE2    100    -    ENSG2    ENST3
```

## `perover`

BED selector input:

```text
chr1    100    180    active1
chr1    240    270    active2
```

Headered or headerless text selector input:

```text
region
chr1:100-180
chr1/240/270
```

```bash
sjcab-peak2anno-db filter-bed -b genes.bed -m perover -i active.txt -o .
```

Result:

```text
chr1    100    300    GENE1    180    +    ENSG1    ENST1
```

## `isoID`

The selector is a transcript ID list:

```text
ENST2
ENST3
```

```bash
sjcab-peak2anno-db dedup-bed hg38 v31 -m isoID -b genes.bed -i ids.txt -o .
```

Result:

```text
chr1    100    280    GENE1    160    +    ENSG1    ENST2
chr1    500    620    GENE2    100    -    ENSG2    ENST3
```

## `isoexp`

The selector is a two-column transcript ID and expression table:

```text
transcript_id    expression
ENST1            2.0
ENST2            8.5
ENST3            1.0
```

```bash
sjcab-peak2anno-db dedup-bed hg38 v31 -m isoexp -b genes.bed -i expression.tsv -o .
```

Result:

```text
chr1    100    280    GENE1    160    +    ENSG1    ENST2
chr1    500    620    GENE2    100    -    ENSG2    ENST3
```

Region text accepts the same common delimiters as `sjcab_peak2anno`, including
`:`, `-`, `_`, `/`, `;`, and `,`.
