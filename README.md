# sjcab_peak2anno_db

Gene-backed peak-to-annotation BED resources for SJ/CAB workflows.

The package installs bundled resources and can also download or regenerate
GENCODE, Roadmap ChromHMM, blacklist, and CpG island BED files. Generated user
data defaults to `~/.sjcab_peak2anno_db`; override it with
`SJCAB_PEAK2ANNO_DB_PATH` or with command-level `--data-dir`/`-d` flags.

```bash
pip install sjcab_peak2anno_db

export SJCAB_PEAK2ANNO_DB_PATH=/path/to/sjcab_peak2anno_db
sjcab-peak2anno-db install
sjcab-peak2anno-db install gencode-bed gencode-feature
sjcab-peak2anno-db install --overwrite
sjcab-peak2anno-db list
```

Supported bundled gene species are `hg19`, `hg38`, `mm10`, `mm9`, `mm39`, and
`sacCer3`. CGI resources are available for `hg38`, `hg19`, `mm10`, `mm9`, and
`mm39`.

Every command that downloads a file appends the source URL to
`{data_dir}/download_urls.log`. Long-running install/download commands print
progress to stderr; stdout is reserved for generated paths and command results.

## GENCODE BED Related

GENCODE BED resources provide transcript-level gene annotations. Two isoform
sets are generated:

- `all`: all transcript isoforms.
- `deduplong`: one longest isoform per gene, grouped by gene symbol by default.

The `def` directory points to the registry-defined default version for each
species, not necessarily the newest parsed GENCODE release.

### `install-bed` / `download-bed`

`install-gencode-bed` installs the package's configured GENCODE BED set into the
user data directory. `download-gencode-bed` processes one species/version into a
chosen output directory.

```bash
sjcab-peak2anno-db install-bed
sjcab-peak2anno-db install-bed -d /path/to/sjcab_peak2anno_db
sjcab-peak2anno-db install-bed --overwrite

sjcab-peak2anno-db download-bed hg38 v31 -o annotations
sjcab-peak2anno-db download-bed hg19 v31lift37 -o annotations
sjcab-peak2anno-db download-bed mm10 vM22 -g gencode.vM22.annotation.gtf.gz -o annotations
sjcab-peak2anno-db download-bed human 100 --source ensembl -o annotations
sjcab-peak2anno-db download-bed human def --source ensembl -o annotations
```

`install-bed` reuses existing generated files by default. Use
`--overwrite` when you intentionally want to regenerate the whole installed BED
layout.

`download-bed` downloads or reuses the expected GTF. It first checks for
the GTF under the output directory, then under the current working directory,
before downloading. Use `--gtf-path`/`-g` to force a specific local GTF or
`--url`/`-u` to use a custom URL. Use `--source ensembl` for species not covered
by the bundled GENCODE registry. Its version is the Ensembl release number
(`100`, for example), or `def` for the latest GTF listed by Ensembl. The same
option is available on `download-feature`.

Default layout:

```text
{data_dir}/bed/{species}/{version}/all.gene.bed
{data_dir}/bed/{species}/{version}/deduplong.gene.bed
{data_dir}/bed/{species}/def -> {version}
```

TSS and TES BED files are not stored by install/download layouts; derive them
from the selected gene BED with `write_tss` or `write_tes` when needed.

Python API:

```python
import sjcab_peak2anno_db as db

db.download_and_convert_gencode_gtf("hg38", "v31", output_dir="annotations")
db.write_tss("all.gene.bed", "all.tss.bed")
db.write_tes("all.gene.bed", "all.tes.bed")
db.write_deduplong("all.gene.bed", "deduplong.gene.bed")
```

### `dedup-bed` / `filter-bed`

Both commands select one isoform per gene from an all-isoform GENCODE BED and
write `{prefix}.gene.bed`, `{prefix}.tss.bed`, and `{prefix}.tes.bed`.

The default `dedup-gencode-bed` selector is `longcol5`: it selects the isoform
with the largest numeric BED column 5. Use `long` to select by interval length
(`end - start`):

`dedup-gencode-bed` falls back to the longest isoform for genes without selector
support. `filter-gencode-bed` uses the same selector logic but omits genes
without selector support.

```bash
sjcab-peak2anno-db dedup-bed hg38 v31 -b all.gene.bed -o annotations
sjcab-peak2anno-db dedup-bed hg38 v31 -m long -b all.gene.bed -o annotations
sjcab-peak2anno-db dedup-bed hg38 v31 -m peak -i h3k4me3_peaks.bed -o annotations
sjcab-peak2anno-db dedup-bed mm10 vM22 -m isoID -i isoforms.txt -K ensid -o annotations
sjcab-peak2anno-db filter-bed -b annotations/bed/hg38/v31/all.gene.bed -m perover -i active_chromhmm.bed -o annotations
sjcab-peak2anno-db filter-bed hg38 v31 -m isoexp -i isoform_expression.tsv --exclusive -o annotations
```

Selection methods:

- `longcol5`: no selector is needed; select the isoform with the largest
  numeric value in BED column 5. This is the default for `dedup-gencode-bed`.
- `long`: no selector is needed; select the isoform with the largest `end - start`.
- `peak`: selector is a peak BED file with peak score in column 5; the isoform
  whose TSS +/- promoter window has the highest peak score is selected. Text
  selectors are also accepted, for example `chr1:100-200 10`.
- `isoID`: selector is a transcript ID list, one ID per line.
- `isoexp`: selector is a two-column table: transcript ID, then expression.
- `perover`: selector is a BED file, often user-filtered ChromHMM active states;
  the isoform promoter with the highest percent overlap is selected. Text
  selectors are also accepted, for example `chr1:100-200` or `chr1_100_200`.

Defaults and naming:

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

Python API:

```python
import sjcab_peak2anno_db as db

db.dedup_gencode_bed(
    "peak",
    "h3k4me3_peaks.bed",
    output_dir="annotations",
    species="hg38",
    version="v31",
)

db.filter_gencode_bed(
    "perover",
    "active_chromhmm.bed",
    output_dir="annotations",
    gene_bed="annotations/bed/hg38/v31/all.gene.bed",
    promoter_bp="2kb",
)
```

### `path`

Use `path` to print an installed resource path. Add `--install`/`-I` when the
resource should be installed first.

```bash
sjcab-peak2anno-db path hg38 gene
```

Python API:

```python
import sjcab_peak2anno_db as db

print(db.path("hg38", "gene"))
```

## GENCODE Feature

GENCODE feature resources are CAB-style merged region BEDs derived from GTF and
gene BED input. The default prefix is the promoter size label, usually `2kb`.

### `install-feature` / `download-feature`

`install-gencode-feature` installs default feature builds into the user data
directory. `download-gencode-feature` writes feature files into a staging/output
directory.

```bash
sjcab-peak2anno-db install-feature all all
sjcab-peak2anno-db install-feature hg38 v31
sjcab-peak2anno-db install-feature hg38 v31 -o feature_downloads

sjcab-peak2anno-db download-feature hg38 v31 -o feature_downloads
sjcab-peak2anno-db download-feature all all -o feature_downloads
sjcab-peak2anno-db download-feature hg38 v31 -o feature_downloads -p 2kb -D 50kb -e 2kb
```

`install-gencode-feature -o DIR` reuses preprocessed feature files from `DIR`
when `order.lst` and the expected feature BED files are already present. Like
`download-gencode-bed`, feature generation reuses an existing expected GTF from
the output directory or current working directory before downloading.
`install-gencode-feature` reuses existing generated files by default; pass
`--overwrite` to rebuild the installed GENCODE BED prerequisite and feature
files.

Default installed layout:

```text
{data_dir}/feature/{species}/{version}/{prefix}/{prefix}.promoter.up.bed
{data_dir}/feature/{species}/{version}/{prefix}/{prefix}.promoter.down.bed
{data_dir}/feature/{species}/{version}/{prefix}/{prefix}.exon.bed
{data_dir}/feature/{species}/{version}/{prefix}/{prefix}.intron.bed
{data_dir}/feature/{species}/{version}/{prefix}/{prefix}.tes.bed
{data_dir}/feature/{species}/{version}/{prefix}/{prefix}.dis5.bed
{data_dir}/feature/{species}/{version}/{prefix}/{prefix}.dis3.bed
{data_dir}/feature/{species}/{version}/{prefix}/{prefix}.intergenic.bed
{data_dir}/feature/{species}/{version}/{prefix}/order.lst
{data_dir}/feature/{species}/def -> {version}/{prefix}
```

`order.lst` contains:

```text
promoter.up
promoter.down
exon
intron
tes
dis5
dis3
intergenic
```

Main feature parameters:

- `--promoter-bp`/`-p`: promoter flank size, default `2kb`.
- `--distal-bp`/`-D`: distal flank size, default `50kb`.
- `--tes-bp`/`-e`: TES flank size, default `2kb`.
- `--prefix`/`-P`: output prefix, default is the promoter size label.
- `--gene-bed`/`-b`: use an existing gene BED for region generation.
- `--gtf-path`/`-g`: use an existing local GTF.

Python API:

```python
import sjcab_peak2anno_db as db

db.download_gencode_feature("hg38", "v31", "feature_downloads")

db.write_tss_flank_region_unions(
    "annotations/bed/hg38/v31/all.gene.bed",
    "feature_downloads",
    promoter_bp="2kb",
    distal_bp="50kb",
)
```

## ChromHMM

ChromHMM commands download Roadmap dense state BED files and write normalized
metadata.

### `install-chromhmm` / `download-chromhmm`

`install-chromhmm` writes into the user data directory. `download-chromhmm`
writes into the selected output directory.

```bash
sjcab-peak2anno-db install-chromhmm -m 18 -G hg19 -i E001,E063
sjcab-peak2anno-db install-chromhmm -m 18 -G hg38 -i E063

sjcab-peak2anno-db download-chromhmm -o chromhmm_downloads -m 25 -G hg19 -c GM12878
sjcab-peak2anno-db download-chromhmm -o chromhmm_downloads -m 18 -G hg38 -i E063
sjcab-peak2anno-db download-chromhmm -o chromhmm_downloads -m 15 -t brain
```

ChromHMM defaults to Roadmap hg19 BED files. Use `--genome hg38`/`-G hg38` for
Roadmap lifted-over BEDs. Select records with epigenome IDs (`--ids`/`-i`),
fuzzy tissue/group text (`--tissue`/`-t`), or fuzzy sample/cell-line text
(`--cellline`/`-c`).

Output layout:

```text
{data_dir}/chromhmm/{genome}/{model}state/{epigenome_id}_{model}_{model_label}_dense.bed.gz
{data_dir}/chromhmm/metadata.tsv
```

`metadata.tsv` includes the downloaded path, genome (`hg19` or `hg38`), state
model, Roadmap epigenome ID, and normalized sample metadata.

Python API:

```python
import sjcab_peak2anno_db as db

db.download_chromhmm(model=18, genome="hg38", ids="E063")
db.download_chromhmm(model=25, genome="hg19", cellline="GM12878")
```

## Segway

Segway commands download the hg19 Segway encyclopedia of human regulatory
elements by trying the ENCODE publication page first
(`https://www.encodeproject.org/publications/94941f71-80c8-43d2-809b-25161efc3be0/`),
then falling back to `https://noble.gs.washington.edu/proj/encyclopedia/`.
If ENCODE returns a human-verification page or no Segway download links, the CLI
prints a fallback message and records only the source that actually provided the
downloaded files.

### `install-segway` / `download-segway`

`install-segway` writes into the user data directory. `download-segway` writes
into the selected output directory. With no selector, the command downloads
`segway_encyclopedia.bed.gz`. Existing files are skipped by default; pass
`--overwrite` to replace them.

```bash
sjcab-peak2anno-db install-segway
sjcab-peak2anno-db install-segway -i GM12878,H1-HESC --include-label-info

sjcab-peak2anno-db download-segway -o segway_downloads -i GM12878
sjcab-peak2anno-db download-segway -o segway_downloads -t brain
sjcab-peak2anno-db download-segway -o segway_downloads --include-caas
```

Segway source files are hg19 only. If another genome is requested, the CLI asks
whether to download hg19 and write a CrossMap liftover helper script. In
non-interactive runs, use `--yes-liftover`.

```bash
sjcab-peak2anno-db download-segway -o segway_downloads -G hg38 --yes-liftover -i GM12878
```

The generated script first tries to activate an existing `segway-liftover`
environment with micromamba, conda, or mamba. If activation succeeds, the script
does not install or modify packages inside it. If activation fails for all
available tools, the script creates the environment with CrossMap. It downloads
the UCSC `hg19To<Genome>.over.chain.gz` file from
`https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/` and lifts over
downloaded `.bed.gz` files. For `install-segway`, lifted files are staged and
then copied into the installed Segway cache under
`~/.sjcab_peak2anno_db/segway/{genome}` unless `-d` selects another data
directory.

Output layout:

```text
{data_dir}/segway/hg19/segway_encyclopedia.bed.gz
{data_dir}/segway/hg19/caas.bed.gz
{data_dir}/segway/hg19/label_info.tab
{data_dir}/segway/hg19/{sample}.bed.gz
{data_dir}/segway/metadata.tsv
{data_dir}/segway/liftover_hg19_to_{genome}.sh
```

Lifted BED files include `lift` in the genome suffix, for example:

```text
{data_dir}/segway/hg38/PERIPHERAL_BLOOD_MONONUCLEAR_PRIMARY_CELLS.hg38lift.bed.gz
```

Python API:

```python
import sjcab_peak2anno_db as db

db.download_segway(names="GM12878,H1-HESC")
db.download_segway(tissue="brain")
db.write_segway_liftover_script(target_genome="hg38")
```

## Blacklists And CGI

Blacklist files are bundled. CGI files are bundled for install and can also be
freshly downloaded from UCSC `cpgIslandExt`.

### `install-blacklists`

`install-blacklists` installs bundled blacklist BED files into the user data
directory. There is no separate blacklist download command.

```bash
sjcab-peak2anno-db install-blacklists
sjcab-peak2anno-db install blacklists
```

Output layout:

```text
{data_dir}/blacklists/{name}.bed.20230411
{data_dir}/blacklists/{name}.bed
```

The current `{name}.bed` path is refreshed as a symlink when supported by the
filesystem.

Python API:

```python
import sjcab_peak2anno_db as db

db.install_blacklists()
```

### `install-cgi` / `download-cgi`

`install-cgi` installs packaged CGI BED files. `download-cgi` downloads fresh
UCSC CpG island tables and writes BED files.

```bash
sjcab-peak2anno-db install-cgi
sjcab-peak2anno-db install cgi

sjcab-peak2anno-db download-cgi
sjcab-peak2anno-db download-cgi -d cgi_downloads --species hg38
```

Output layout:

```text
{data_dir}/cgi/{species}_cgi.bed
```

Python API:

```python
import sjcab_peak2anno_db as db

db.install_cgi()
db.download_cgi(species="hg38")
```

## Command Examples

```bash
sjcab-peak2anno-db list
sjcab-peak2anno-db install
sjcab-peak2anno-db install gencode-bed gencode-feature blacklists cgi
sjcab-peak2anno-db install -c gencode-bed -c gencode-feature
sjcab-peak2anno-db install --overwrite

sjcab-peak2anno-db install-bed
sjcab-peak2anno-db download-bed hg38 v31 -o annotations
sjcab-peak2anno-db download-bed mm10 vM22 -g gencode.vM22.annotation.gtf.gz -o annotations
sjcab-peak2anno-db dedup-bed hg38 v31 -m peak -i h3k4me3_peaks.bed -o annotations
sjcab-peak2anno-db dedup-bed mm10 vM22 -m isoID -i isoforms.txt -K ensid -o annotations
sjcab-peak2anno-db filter-bed -b annotations/bed/hg38/v31/all.gene.bed -m perover -i active_chromhmm.bed -o annotations
sjcab-peak2anno-db filter-bed hg38 v31 -m isoexp -i isoform_expression.tsv --exclusive -o annotations

sjcab-peak2anno-db install-feature all all
sjcab-peak2anno-db install-feature hg38 v31 -o feature_downloads
sjcab-peak2anno-db download-feature hg38 v31 -o feature_downloads
sjcab-peak2anno-db download-feature all all -o feature_downloads

sjcab-peak2anno-db install-chromhmm -m 18 -G hg19 -i E001,E063
sjcab-peak2anno-db install-chromhmm -m 18 -G hg38 -i E063
sjcab-peak2anno-db download-chromhmm -o chromhmm_downloads -m 25 -G hg19 -c GM12878
sjcab-peak2anno-db download-chromhmm -o chromhmm_downloads -m 18 -G hg38 -i E063

sjcab-peak2anno-db install-segway -i GM12878,H1-HESC
sjcab-peak2anno-db download-segway -o segway_downloads -t brain
sjcab-peak2anno-db download-segway -o segway_downloads -G hg38 --yes-liftover -i GM12878

sjcab-peak2anno-db install-blacklists
sjcab-peak2anno-db install-cgi
sjcab-peak2anno-db download-cgi
sjcab-peak2anno-db download-cgi -d cgi_downloads --species hg38

sjcab-peak2anno-db path hg38 gene
```
