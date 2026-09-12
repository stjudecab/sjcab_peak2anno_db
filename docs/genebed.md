# Genebed resources

Gene BED resources are transcript-level annotations. `all` contains every
transcript isoform; `deduplong` keeps one longest isoform per gene. The
`def` link points to the registry-defined default version.

Bundled GENCODE builds include `hg19`, `hg38`, `mm9`, `mm10`, and `mm39`.
Other assemblies are resolved from GENCODE, UCSC, or Ensembl.

## Commands

```bash
sjcab-peak2anno-db install-genebed
sjcab-peak2anno-db install-genebed mm10 vM22 -g gencode.vM22.annotation.gtf.gz
sjcab-peak2anno-db download-genebed hg38 v31 -o annotations
sjcab-peak2anno-db download-genebed hg19 v31lift37 -o annotations
```

For Ensembl, `def`, `default`, `current`, and `latest` resolve the current
release. Vertebrate and Ensembl Genomes catalogues are cached by release;
`def` links to the current release. Assembly, species, name, and
`assembly_accession` can be used as species identifiers when unique.

For UCSC builds, GTF metadata and liftOver chain URLs are cached in
`{db-path}/ucsc/gtf_builds.tsv`. Sizes files use the assembly name whenever
possible; `.sizes.clean` uses primary chromosomes from the assembly report.

## Layout

```text
{db-path}/bed/{species}/{version}/all.gene.bed
{db-path}/bed/{species}/{version}/deduplong.gene.bed
{db-path}/bed/{species}/def -> {version}
{db-path}/sizes/{species}.sizes
{db-path}/sizes/{species}.sizes.clean
```

## Feature generation from a local GTF

```bash
sjcab-peak2anno-db download-genebed mm10 vM22 \
  --gtf-path gencode.vM22.annotation.gtf.gz \
  --output-dir annotations \
  --processes 8
```

Use `--gene-type` repeatedly to restrict GENCODE types. See
[Feature BEDs](features.md) for derived regions and
[Deduplication](deduplication.md) for transcript selection.
