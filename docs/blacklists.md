# Blacklists

Install the bundled blacklist BED files:

```bash
sjcab-peak2anno-db install-blacklists
sjcab-peak2anno-db install blacklists
```

Outputs are stored as `{db-path}/blacklists/{name}.bed`. For provenance and
source details, see the [bundled blacklist notes](https://github.com/stjudecab/sjcab_peak2anno_db/blob/main/src/sjcab_peak2anno_db/data/blacklists/readme.md).

## Resource URLs

- [Boyle Lab blacklist lists](https://github.com/Boyle-Lab/Blacklist/tree/master/lists)
- [ENCODE blacklist annotation](https://www.encodeproject.org/annotations/ENCSR636HFF/)
- [Kundaje ENCODE blacklist notes](https://personal.broadinstitute.org/anshul/projects/encode/rawdata/blacklists/hg19-blacklist-README.pdf)

## Liftover comparison

The following is an interval-overlap comparison of the generated liftover
files against the native target blacklist files:

- `TP`: native truth intervals with at least one base overlapping a liftover
  interval.
- `FN`: native truth intervals without an overlap (`native truth - TP`).
- `TPR = TP / native truth`; `FNR = FN / native truth`.
- `merge` is applied to the liftover result before overlap testing. `none`
  means sorted but unmerged input; the other values are `bedtools merge -d N`.

| assembly | liftover | merge | native truth | liftover intervals | TP | TPR | FN | FNR |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| hg19 | CrossMap | none | 1,871 | 59,857 | 1,426 | 76.216% | 445 | 23.784% |
| hg19 | CrossMap | 100 | 1,871 | 7,713 | 1,426 | 76.216% | 445 | 23.784% |
| hg19 | CrossMap | 200 | 1,871 | 6,575 | 1,426 | 76.216% | 445 | 23.784% |
| hg19 | UCSC liftOver | none | 1,871 | 2,237 | 1,255 | 67.076% | 616 | 32.924% |
| hg19 | UCSC liftOver | 100 | 1,871 | 2,186 | 1,255 | 67.076% | 616 | 32.924% |
| hg19 | UCSC liftOver | 200 | 1,871 | 2,158 | 1,255 | 67.076% | 616 | 32.924% |
| mm10 | CrossMap | none | 4,874 | 270,480 | 64 | 1.313% | 4,810 | 98.687% |
| mm10 | CrossMap | 100 | 4,874 | 5,478 | 64 | 1.313% | 4,810 | 98.687% |
| mm10 | CrossMap | 200 | 4,874 | 4,828 | 64 | 1.313% | 4,810 | 98.687% |
| mm10 | UCSC liftOver | none | 4,874 | 232 | 3 | 0.062% | 4,871 | 99.938% |
| mm10 | UCSC liftOver | 100 | 4,874 | 69 | 3 | 0.062% | 4,871 | 99.938% |
| mm10 | UCSC liftOver | 200 | 4,874 | 68 | 3 | 0.062% | 4,871 | 99.938% |

### Prediction-side overlap

If the question is “how many liftover regions overlap native truth?”, the
corresponding metric is prediction precision/positive predictive value (PPV),
not TPR. Here `liftover TP` is a liftover interval overlapping at least one
native truth interval, and `liftover FP` is a liftover interval with no native
overlap. The native-truth recall columns are retained for reference.

| assembly | liftover | merge | liftover intervals | liftover TP | PPV | liftover FP | native TP | native TPR |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| hg19 | CrossMap | none | 59,857 | 47,262 | 78.958% | 12,595 | 1,426 | 76.216% |
| hg19 | CrossMap | 100 | 7,713 | 3,286 | 42.603% | 4,427 | 1,426 | 76.216% |
| hg19 | CrossMap | 200 | 6,575 | 2,773 | 42.175% | 3,802 | 1,426 | 76.216% |
| hg19 | UCSC liftOver | none | 2,237 | 1,342 | 59.991% | 895 | 1,255 | 67.076% |
| hg19 | UCSC liftOver | 100 | 2,186 | 1,329 | 60.796% | 857 | 1,255 | 67.076% |
| hg19 | UCSC liftOver | 200 | 2,158 | 1,329 | 61.585% | 829 | 1,255 | 67.076% |
| mm10 | CrossMap | none | 270,480 | 4,487 | 1.659% | 265,993 | 64 | 1.313% |
| mm10 | CrossMap | 100 | 5,478 | 150 | 2.738% | 5,328 | 64 | 1.313% |
| mm10 | CrossMap | 200 | 4,828 | 134 | 2.776% | 4,694 | 64 | 1.313% |
| mm10 | UCSC liftOver | none | 232 | 4 | 1.724% | 228 | 3 | 0.062% |
| mm10 | UCSC liftOver | 100 | 69 | 4 | 5.797% | 65 | 3 | 0.062% |
| mm10 | UCSC liftOver | 200 | 68 | 4 | 5.882% | 64 | 3 | 0.062% |

### Confidence and limitations

These figures are overlap-based recall only. They do not measure false
positive intervals, coordinate accuracy, one-to-many mappings, or whether the
liftover output retains the intended blacklist boundaries. The very large
CrossMap interval counts, particularly for mm10, indicate fragmentation or
one-to-many mappings and should not be interpreted as higher confidence. The
mm10 blacklist results have extremely low recovery with both methods and
should be treated as unsuitable without additional review against a native
mm10 blacklist release. Results can also change with the selected chain file,
input release, chromosome naming, and merge policy.
