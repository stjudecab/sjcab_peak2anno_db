# CGI

Install packaged CpG-island BEDs or download fresh UCSC `cpgIslandExt`
tables:

```bash
sjcab-peak2anno-db install-cgi
sjcab-peak2anno-db install cgi
sjcab-peak2anno-db download-cgi --species hg38 -d cgi_downloads
```

Outputs are stored as `{db-path}/cgi/{species}_cgi.bed`.

## Resource URLs

- [UCSC genome downloads](https://hgdownload.soe.ucsc.edu/downloads.html)
- [UCSC per-assembly database tables](https://hgdownload.soe.ucsc.edu/goldenPath/)
- [UCSC CpG island table pattern](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz)

## Liftover comparison

The following is an interval-overlap comparison of the generated liftover
files against native target CGI BED files:

- The native intervals are the true/reference set; liftover intervals are the
  predictions being evaluated.
- `TP`: native truth intervals with at least one base overlapping a liftover
  interval.
- `FN`: native truth intervals without an overlap (`native truth - TP`).
- `TPR = TP / native truth`; `FNR = FN / native truth`.
- `merge` is applied to the liftover result before overlap testing. `none`
  means sorted but unmerged input; the other values are `bedtools merge -d N`.

| assembly | liftover | merge | native truth | liftover intervals | TP | TPR | FN | FNR |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| hg19 | CrossMap | none | 30,344 | 33,643 | 27,183 | 89.583% | 3,161 | 10.417% |
| hg19 | CrossMap | 100 | 30,344 | 27,777 | 27,183 | 89.583% | 3,161 | 10.417% |
| hg19 | CrossMap | 200 | 30,344 | 27,555 | 27,183 | 89.583% | 3,161 | 10.417% |
| hg19 | UCSC liftOver | none | 30,344 | 29,997 | 27,112 | 89.349% | 3,232 | 10.651% |
| hg19 | UCSC liftOver | 100 | 30,344 | 27,605 | 27,112 | 89.349% | 3,232 | 10.651% |
| hg19 | UCSC liftOver | 200 | 30,344 | 27,411 | 27,112 | 89.349% | 3,232 | 10.651% |
| mm10 | CrossMap | none | 17,017 | 479,438 | 13,735 | 80.713% | 3,282 | 19.287% |
| mm10 | CrossMap | 100 | 17,017 | 23,535 | 13,735 | 80.713% | 3,282 | 19.287% |
| mm10 | CrossMap | 200 | 17,017 | 22,755 | 13,735 | 80.713% | 3,282 | 19.287% |
| mm10 | UCSC liftOver | none | 17,017 | 5,272 | 2,523 | 14.826% | 14,494 | 85.174% |
| mm10 | UCSC liftOver | 100 | 17,017 | 4,717 | 2,523 | 14.826% | 14,494 | 85.174% |
| mm10 | UCSC liftOver | 200 | 17,017 | 4,698 | 2,523 | 14.826% | 14,494 | 85.174% |

### Prediction-side overlap

If the question is “how many liftover regions overlap native truth?”, the
corresponding metric is prediction precision/positive predictive value (PPV),
not TPR. Here `liftover TP` is a liftover interval overlapping at least one
native truth interval, and `liftover FP` is a liftover interval with no native
overlap. The native-truth recall columns are retained for reference.

| assembly | liftover | merge | liftover intervals | liftover TP | PPV | liftover FP | native TP | native TPR |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| hg19 | CrossMap | none | 33,643 | 32,329 | 96.094% | 1,314 | 27,183 | 89.583% |
| hg19 | CrossMap | 100 | 27,777 | 27,037 | 97.336% | 740 | 27,183 | 89.583% |
| hg19 | CrossMap | 200 | 27,555 | 26,835 | 97.387% | 720 | 27,183 | 89.583% |
| hg19 | UCSC liftOver | none | 29,997 | 29,283 | 97.620% | 714 | 27,112 | 89.349% |
| hg19 | UCSC liftOver | 100 | 27,605 | 26,934 | 97.569% | 671 | 27,112 | 89.349% |
| hg19 | UCSC liftOver | 200 | 27,411 | 26,748 | 97.581% | 663 | 27,112 | 89.349% |
| mm10 | CrossMap | none | 479,438 | 251,575 | 52.473% | 227,863 | 13,735 | 80.713% |
| mm10 | CrossMap | 100 | 23,535 | 13,392 | 56.903% | 10,143 | 13,735 | 80.713% |
| mm10 | CrossMap | 200 | 22,755 | 13,222 | 58.106% | 9,533 | 13,735 | 80.713% |
| mm10 | UCSC liftOver | none | 5,272 | 2,734 | 51.859% | 2,538 | 2,523 | 14.826% |
| mm10 | UCSC liftOver | 100 | 4,717 | 2,414 | 51.177% | 2,303 | 2,523 | 14.826% |
| mm10 | UCSC liftOver | 200 | 4,698 | 2,409 | 51.277% | 2,289 | 2,523 | 14.826% |

### Confidence and limitations

These figures are overlap-based recall only. They do not measure false
positive intervals, coordinate accuracy, one-to-many mappings, or preservation
of CGI metadata. The very large CrossMap interval counts, especially for
mm10, indicate fragmentation or one-to-many mappings and should not be
interpreted as higher confidence. UCSC liftOver recovers only 14.826% of the
mm10 native intervals in this comparison, so that result should not be used
without review against a native mm10 CGI release. Results can change with the
selected chain file, source release, chromosome naming, and merge policy.
