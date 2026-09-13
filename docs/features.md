# FeatureBEDs

Feature resources are merged CAB-style regions derived from a GTF and gene
BED. The default output prefix is the promoter size, normally `2kb`.

## Commands

```bash
sjcab-peak2anno-db install-feature hg38 v31
sjcab-peak2anno-db download-feature hg38 v31 -o feature_downloads
sjcab-peak2anno-db download-feature hg38 v31 \
  -p 2kb -D 50kb -e 2kb -j 8
```

Omitting the version for `install-feature` uses `def`. `-name` stores a species
under a custom name and records the assembly-to-name mapping in
`custom.name.tsv`. `-dry-run` tests and prints the final non-404 GTF URL.

Main options:

- `--promoter-bp`/`-p`: promoter flank, default `2kb`.
- `--distal-bp`/`-D`: distal flank, default `50kb`.
- `--tes-bp`/`-e`: TES flank, default `2kb`.
- `--prefix`/`-P`: output prefix.
- `--gene-bed`/`-b`: existing GeneBED.
- `--gtf-path`/`-g`: existing local GTF.
- `--processes`/`-j`: parallel conversion and feature-file writing.

## Output layout

```text
{db-path}/feature/{species}/{version}/{prefix}/{prefix}.promoter.up.bed
{db-path}/feature/{species}/{version}/{prefix}/{prefix}.promoter.down.bed
{db-path}/feature/{species}/{version}/{prefix}/{prefix}.exon.bed
{db-path}/feature/{species}/{version}/{prefix}/{prefix}.intron.bed
{db-path}/feature/{species}/{version}/{prefix}/{prefix}.tes.bed
{db-path}/feature/{species}/{version}/{prefix}/{prefix}.dis5.bed
{db-path}/feature/{species}/{version}/{prefix}/{prefix}.dis3.bed
{db-path}/feature/{species}/{version}/{prefix}/{prefix}.intergenic.bed
{db-path}/feature/{species}/{version}/{prefix}/order.lst
{db-path}/feature/{species}/def -> {version}/{prefix}
```

`order.lst` contains `promoter.up`, `promoter.down`, `exon`, `intron`, `tes`,
`dis5`, `dis3`, and `intergenic` in that order.
