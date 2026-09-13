# Deduplication and selectors

`dedup-bed` selects one transcript per gene and falls back to the longest
transcript when no selector matches. `filter-bed` uses the same selectors but
omits genes without a selector match.

```bash
sjcab-peak2anno-db dedup-bed hg38 v31 -m longcol5 -b all.gene.bed -o annotations
sjcab-peak2anno-db dedup-bed hg38 v31 -m long -b all.gene.bed -o annotations
sjcab-peak2anno-db dedup-bed hg38 v31 -m peak -i h3k4me3_peaks.bed -o annotations
sjcab-peak2anno-db dedup-bed mm10 vM22 -m isoID -i isoforms.txt -K ensid -o annotations
sjcab-peak2anno-db filter-bed -m perover -i active_chromhmm.bed -o annotations
```

## Selector methods

- `longcol5`: largest numeric BED column 5; the default for `dedup-bed`.
- `long`: largest genomic interval, `end - start`.
- `peak`: highest peak score in the transcript promoter window.
- `isoID`: transcript ID list, one ID per line.
- `isoexp`: two-column transcript ID and expression table.
- `perover`: greatest overlap with a selector BED file.

The promoter half-window defaults to `2kb` and can be changed with
`--promoter-bp`/`-p`. Use `--exclusive` for exact transcript IDs or fully
contained BED features.

Use `--promoter-down` to set a separate downstream promoter window; it defaults
to `--promoter-bp`. Selector files are auto-detected once per file from the
first data row. BED selectors use score column 5 and text selectors use score
column 2 by default. Configure the columns with
`SJCAB_PEAK2ANNO_DB_BED_SCORE_COLUMN` and
`SJCAB_PEAK2ANNO_DB_TXT_SCORE_COLUMN`. Configure delimiters inside text-mode
peak regions with `SJCAB_PEAK2ANNO_DB_PEAK_TXT_DELIMITER`.

By default transcripts are grouped by gene symbol. Use `--gene-key ensid` to
group by Ensembl/GENCODE gene ID.

## Python API

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
)
```

For complete selector input/output examples, see the repository's
[README.DEDUP.md](https://github.com/stjudecab/sjcab_peak2anno_db/blob/main/README.DEDUP.md).
