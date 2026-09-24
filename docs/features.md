# FeatureBEDs

Feature resources are merged CAB-style regions derived from a GTF and gene
BED. The default output prefix is the promoter size, normally `2kb`.

## Commands

```bash
sjcab-peak2anno-db install-feature hg38 v31
sjcab-peak2anno-db download-feature hg38 v31 -o feature_downloads
sjcab-peak2anno-db download-feature hg38 v31 \
  -p 2k -D 50kb -e 2k
```

Omitting the version for `install-feature` uses `def`. `-name` stores a species
under a custom name and records the assembly-to-name mapping in
`custom.name.tsv`. `-dry-run` tests and prints the final non-404 GTF URL.

Main options:

- `--promoter-bp`/`-p`: promoter flank, default `2kb`.
- `--promoter-down`: downstream promoter flank; defaults to `--promoter-bp`.
- `--distal-bp`/`-D`: distal flank, default `50kb`.
- `--distal-down`: downstream distal flank; defaults to `--distal-bp`.
- `--tes-bp`/`-e`: TES flank, default `2kb`.
- `--tes-up`: upstream TES flank; defaults to `--tes-bp`.
- `--prefix`/`-P`: output prefix.
- `--gene-bed`/`-b`: existing GeneBED.
- `--gtf-path`/`-g`: existing local GTF.
- `--processes`/`-j`: parallel processing of independent GTFs. Downloads are
  sequential and conversion within each worker is single-process. The worker
  count is capped at the CPUs allocated to the job.

## Output layout

```text
{db-path}/feature/{species}/{version}/{prefix}/{prefix}.promoter.up.bed
{db-path}/feature/{species}/{version}/{prefix}/{prefix}.5utr.bed
{db-path}/feature/{species}/{version}/{prefix}/{prefix}.3utr.bed
{db-path}/feature/{species}/{version}/{prefix}/{prefix}.promoter.down.bed
{db-path}/feature/{species}/{version}/{prefix}/{prefix}.exon.bed
{db-path}/feature/{species}/{version}/{prefix}/{prefix}.intron.bed
{db-path}/feature/{species}/{version}/{prefix}/{prefix}.tes.bed
{db-path}/feature/{species}/{version}/{prefix}/{prefix}.dis5.bed
{db-path}/feature/{species}/{version}/{prefix}/{prefix}.dis3.bed
{db-path}/feature/{species}/{version}/{prefix}/{prefix}.intergenic.bed
{db-path}/feature/{species}/{version}/{prefix}/order.lst
{db-path}/feature/{species}/{version}/{prefix}/order.utr.lst
{db-path}/feature/{species}/def -> {version}/{prefix}
```

`order.lst` and `order.utr.lst` contain tab-separated rows with the BED
filename, feature name, and full name. Both use this order: `Promoter.Up`,
`5UTR`, `3UTR`, `Promoter.Down`, `Exon`, `Intron`, `TES`, `Dis5`, `Dis3`, and
`Intergenic`. For example, the TES row is:

```text
{prefix}.tes.bed  TES  Transcription_End_Sites
```
