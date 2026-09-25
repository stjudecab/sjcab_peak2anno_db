# Changelog

## Unreleased

## 0.2.0

- Ensure packaged BED resources are sorted by chromosome, start, and end.
- Sort generated GeneBED, FeatureBED, TSS/TES, deduplicated, CGI, and
  blacklist outputs regardless of the selected interval-writing backend.
- Add coverage for sorted interval generation and packaged BED validation.
- Complete existing RC files with missing commented settings and normalize
  commented variables to the `#SJCAB_...` form.
- Add aligned component-aware `list` output with a `source` column.
- Store installed GeneBEDs under `genebed/` instead of `bed/`.
- Record installed FeatureBED source GTF filenames in `installed.tsv` and show
  them in `list feature`.

## 0.1.9

- See the release history in the repository for the preceding configuration,
  documentation, and publishing updates.
