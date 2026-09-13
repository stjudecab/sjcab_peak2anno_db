# sjcab_peak2anno_db

Gene-backed peak-to-annotation BED resources for SJ/CAB workflows.

The package installs bundled resources and can also download/generate
from GENCODE/UCSC/Ensembl gtf, Roadmap ChromHMM, blacklist, and CpG island
BED files. Generated user data defaults to `~/.sjcab_peak2anno_db`,
override it with `SJCAB_PEAK2ANNO_DB_PATH` 
or with command-level `--db-path`/`-d` flags.

```bash
# Install the package and its default resource set.
pip install sjcab_peak2anno_db
conda install stjudecab::sjcab_peak2anno_db

export SJCAB_PEAK2ANNO_DB_PATH=~/.sjcab_peak2anno_db
sjcab-peak2anno-db install # install default for "GeneBEDs FeatureBEDs blacklists cgi"
sjcab-peak2anno-db install genebed feature # install selected default group
sjcab-peak2anno-db list # list default bundle

sjcab-peak2anno-db download-genebed hg38 v31 -o annotations # download Gencode hg38 v31 GeneBEDs to annotation folder
sjcab-peak2anno-db install-genebed hg19 v31lift37 # install Gencode hg19 v31lift37 to --db-path
sjcab-peak2anno-db download-feature mm10 vM22 -g gencode.vM22.annotation.gtf.gz -o annotations # convert own gtf file to GeneBED and FeatureBEDs
sjcab-peak2anno-db install-feature cat 116 # install GeneBEDs/FeatureBEDs for cat ensemble release 116
sjcab-peak2anno-db install-feature equcab2 # install GeneBEDs/FeatureBEDs for horse UCSC version equCab2
sjcab-peak2anno-db install-feature zea_mays # install GeneBEDs/FeatureBEDs for corn ensemble current release(empty/def/default/current/latest are the same)
```

Supported bundled GeneBED species are `hg19`, `hg38`, `mm10`, `mm9`,
`mm39`(Gencode) and `sacCer3`. All other GeneBED/FeatureBED resources for species/version are automatically generated from GENCODE/UCSC/Ensembl GTFs.
CGI: `hg38`, `hg19`, `mm10`, `mm9`, and `mm39`
blacklists: [readme.md](https://github.com/stjudecab/sjcab_peak2anno_db/tree/main/src/sjcab_peak2anno_db/data/blacklists)

All download first checks the db-path/cache folder(> output directory > working directory) before downloading. Every downloaded file appends the source URL to `{db-path}/download_urls.log`

## Command conventions and common options

The general difference between command pairs is the destination:

- `install-*`: install or generate resources under the configured `db_path`
  directory.
- `download-*`: write downloaded/generated resources under `--output-dir` (default workdir).
- `*-feature`: installs/downloads both GeneBEDs and FeatureBEDs
- `blacklists|cgi|chromhmm|segway`: if there isn't available resource for selected
  species from UCSC, liftover script can be provided based on chain file from UCSC.

Common option styles:

- `SPECIES`, `-s`, `--species SPECIES`: select a species or genome build. 
  It would match exact by Gencode > UCSC > Ensemble, if no exact match, will 
  print the matches from Ensembl(assembly > latin name > common name) for you 
  to copy and rerun with.
- `VERSION`, `-v`, `--ver VERSION`: select a release/version, such as `v31`,
  `vM39`, `empty|def|default|current|latest`, where supported. The option form can replace
- `species1:version1,species2:version2`: Species and versions accept comma-separated values 
  or a `.lst`/`.list` file with one choice per row. A two-column list supplies 
  `species version` pairs; `species:version` values are also accepted. 
  A single species or version is broadcast across the other list, while equal-length 
  lists are paired in order.
- `-d`, `--db-path DIR`: override the database/cache directory for commands
  that read or install the database.
- `-name NAME`: store the selected GeneBEDs or FeatureBEDs under a custom
  species name and record the assembly-to-name mapping in `custom.name.tsv`.
- `-o`, `--output-dir DIR`: select a staging/output directory for download
  commands. Default for `download-*` use working directory as default.
- `-n`, `--no-overwrite`: preserve existing files. `--overwrite` rewrites them
  where supported.
- `--clean-cache [DAYS]`: remove cache files older than `DAYS` after BED
  generation. The default is 90 days; using it without a value removes the
  current cached GTF immediately, and negative values have the same immediate
  cleanup behavior.
- `--sizes-clean [0|1]`: create `.sizes.clean` files by default. Use
  `--sizes-clean 0` to disable them. 
- `-j`, `--processes N`: use `N` worker processes for GTF-to-GeneBED
  conversion and FeatureBED generation. Default is `4`.
- `-g`, `--gtf-path FILE`: use an existing local GTF instead of downloading
  one.
- `-u`, `--url URL`: override the resolved GTF URL.
- `--ucsc-source {ens,refseq}`: select the UCSC gene-table source when a UCSC
  build is used.
- `-dry-run`: resolve and print downloadable URLs without downloading, where
  supported.

## User configuration

If no RC file exists, the default XDG RC template is created automatically;
all settings in the generated template are commented out.
RC files are read from `~/.sjcab_peak2anno.rc` and
`$XDG_CONFIG_HOME/sjcab_peak2anno/.sjcab_peak2anno.rc` (the XDG file takes
precedence). Set `SJCAB_PEAK2ANNO_CONFIG` to add a specific RC file; that file
has the highest RC precedence. Explicit command-line arguments take precedence
over configuration; environment variables beginning with
`SJCAB_PEAK2ANNO_DB_` take precedence over RC values.

Example:

```text
SJCAB_PEAK2ANNO_DB_PATH=~/.sjcab_peak2anno_db
SJCAB_PEAK2ANNO_DB_INSTALL_OPTIONS=genebed,feature,blacklists,cgi
SJCAB_PEAK2ANNO_DB_INSTALL_SPECIES=hg38,hg19,mm10,mm39
SJCAB_PEAK2ANNO_DB_INSTALL_VERSIONS=v31,v31lift37,vM22,vM39
SJCAB_PEAK2ANNO_DB_VERSION_STALE_DAYS=90
SJCAB_PEAK2ANNO_DB_SIZESCLEAN=1
SJCAB_PEAK2ANNO_DB_CLEANCACHE=90
SJCAB_PEAK2ANNO_DB_TXT_DELIMITER=\s+
SJCAB_PEAK2ANNO_DB_BED_SCORE_COLUMN=5
SJCAB_PEAK2ANNO_DB_TXT_SCORE_COLUMN=2
```

RC keys use the same names as the environment variables; the
`SJCAB_PEAK2ANNO_DB_` prefix is optional for backwards-compatible short keys.
`SJCAB_PEAK2ANNO_DB_PATH` sets the default data directory. `INSTALL_OPTIONS`
accepts `genebed`, `feature`, `blacklists`, and `cgi`.
`INSTALL_SPECIES` and `INSTALL_VERSIONS` control the default species/version
pairs used by the feature installer. They accept comma-separated values or
`.lst`/`.list` files, with single-value broadcasting and equal-length pairing.

`VERSION_STALE_DAYS` controls refresh of the Ensembl `VERSION`,
`species_EnsemblVertebrates.txt`, and `species.txt` catalogs. `SIZESCLEAN=1`
creates `.sizes.clean` files by default; set it to `0` to disable them. The
`--sizes-clean 0` option is available on feature-generation and `install`
commands. `CLEANCACHE` is an age
in days: the default is `90`, so cache files older than 90 days are removed
after BED generation. A negative value removes the newly used cached GTF
immediately. The command-line form is `--clean-cache [DAYS]`; using
`--clean-cache` without a value means immediate cleanup.

`TXT_DELIMITER` controls splitting of text selector rows for `peak` and
`perover`; its commented default is whitespace (`\s+`). Selector files are
auto-detected as BED or text: BED rows use the configured fifth-column score,
while text rows use the configured second-column score. Override these with
`BED_SCORE_COLUMN` and `TXT_SCORE_COLUMN` in the environment or RC file.

## GeneBEDs

GeneBED resources provide transcript-level gene annotations. Two isoform
sets are generated:

- `all`: all transcript isoforms.
- `deduplong`: one longest isoform per gene, grouped by gene symbol by default.

The `def` directory points to the registry-defined default version for each
species, not necessarily the newest parsed GENCODE release.

### `install-genebed` / `download-genebed`

```bash
# Install or regenerate annotation BED resources.
sjcab-peak2anno-db install-genebed

sjcab-peak2anno-db install-genebed mm10 vM22 -g gencode.vM22.annotation.gtf.gz
sjcab-peak2anno-db install-genebed cat
sjcab-peak2anno-db install-genebed dog 100
sjcab-peak2anno-db install-genebed hg38 -name human
sjcab-peak2anno-db download-genebed hg38 v31 -o annotations
sjcab-peak2anno-db download-genebed hg19 v31lift37 -o annotations
```

For Ensembl releases, `def`, `default`, `current`, and `latest` resolve the
current release metadata. The separate links are
`{db-path}/ensembl/vertebrates/def -> {release}` and
`{db-path}/ensembl/genomes/def -> {release}`. Ensembl Genomes determines its
current release from `https://ftp.ebi.ac.uk/pub/ensemblgenomes/VERSION`.

For a species not built into the resolver, the command checks separate cached
catalogs under `{db-path}/ensembl/vertebrates/{release}/` and
`{db-path}/ensembl/genomes/{release}/`, downloading the missing catalog only when
needed. The `def` directory links to the current release. Each catalog stores
`assembly`, `species`, `division`, `name`, and
`assembly_accession` in that order; `assembly` or `species` or `assembly_accession`
can be passed back as the species ID as long as they are unique.
A release-specific reduced catalog is stored at
`{db-path}/ensembl/vertebrates/{release}/` or
`{db-path}/ensembl/genomes/{release}/`; each catalog's `def` link points to
its own release directory.

For automatic GTF selection, the order is GENCODE, a recognized UCSC build
with a UCSC GTF, Ensembl Vertebrates, then Ensembl Genomes. UCSC build IDs and
liftOver chain file URL are cached in `{db-path}/ucsc/gtf_builds.tsv` after reading the
[UCSC downloads page](https://hgdownload.soe.ucsc.edu/downloads.html) and [liftOver](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/) 
would refresh if older than `VERSION_STALE_DAYS`.

`{species}.sizes` for UCSC were based on `bigZips/{build}.chrom.sizes`, otherwise based on ensembl assembly report.
`{species}.sizes.clean` is based on ensembl assembly report only kept primary chromosomes
All generated liftOver scripts store chain files in `{db-path}/cache/chains/` and
reuse an existing non-empty chain file.

Default layout:

```text
{db-path}/bed/{species}/{version}/all.gene.bed
{db-path}/bed/{species}/{version}/deduplong.gene.bed
{db-path}/bed/{species}/def -> {version}
{db-path}/sizes/{species}.sizes
{db-path}/sizes/{species}.sizes.clean
```

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

`dedup-bed` falls back to the longest isoform for genes without selector
support (each gene from original would have at least one transcript).
`filter-bed` uses the same selector logic but omits genes without selector
 support (could remove a lot genes without supports).

The default `dedup-bed` selector is `longcol5`: it selects the isoform
with the largest numeric BED column 5. Use `long` to select by interval length
(`end - start`):

Detailed selector input examples, including input and output BED content, are in
[README.DEDUP.md](https://github.com/stjudecab/sjcab_peak2anno_db/blob/main/README.DEDUP.md).

```bash
# Select one transcript per gene using the default column-5 selector.
sjcab-peak2anno-db dedup-bed hg38 v31 -b all.gene.bed -o annotations
# Select one transcript per gene using genomic interval length.
sjcab-peak2anno-db dedup-bed hg38 v31 -m long -b all.gene.bed -o annotations
# Select transcripts using peak scores around their TSS.
sjcab-peak2anno-db dedup-bed hg38 v31 -m peak -i h3k4me3_peaks.bed -o annotations
# Select transcript IDs explicitly.
sjcab-peak2anno-db dedup-bed mm10 vM22 -m isoID -i isoforms.txt -K ensid -o annotations
# Keep genes whose promoters overlap the selector file.
sjcab-peak2anno-db filter-bed -b annotations/bed/hg38/v31/all.gene.bed -m perover -i active_chromhmm.bed -o annotations
# Select by expression and require exact selector matches.
sjcab-peak2anno-db filter-bed hg38 v31 -m isoexp -i isoform_expression.tsv --exclusive -o annotations
```

Deduplicate Selector methods:

- `longcol5`: no selector is needed; select the isoform with the largest
  numeric value in BED column 5. This is the default for `dedup-bed`.
- `long`: no selector is needed; select the isoform with the largest `end - start`.
- `peak`: selector is a peak BED file with peak score in column 5; the isoform
  whose TSS +/- promoter window has the highest peak score is selected. Text
  selectors are also accepted, for example `chr1:100-200 10`.
- `isoID`: selector is a transcript ID list, one ID per line.
- `isoexp`: selector is a two-column table: transcript ID, then expression.
- `perover`: selector is a BED file, often user-filtered ChromHMM active states;
  the isoform promoter with the highest percent overlap is selected. Text
  selectors are also accepted, for example `chr1:100-200` or `chr1_100_200`.


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

## FeatureBEDs

FeatureBED resources are CAB-style merged region BEDs derived from GTF and GeneBED
input. The default prefix is the promoter size label, usually `2kb`.

### `install-feature` / `download-feature`

```bash
# Install or generate merged FeatureBED resources.
sjcab-peak2anno-db install-feature all all
sjcab-peak2anno-db install-feature hg38 v31
sjcab-peak2anno-db install-feature hg38 v31 -o feature_downloads

sjcab-peak2anno-db download-feature hg38 v31 -o feature_downloads
sjcab-peak2anno-db download-feature all all -o feature_downloads
sjcab-peak2anno-db download-feature hg38 v31 -o feature_downloads -p 2kb -D 50kb -e 2kb
# Resolve and print the GTF URL without downloading.
sjcab-peak2anno-db download-feature horse def -dry-run
# Store one species under a custom name and record the mapping.
sjcab-peak2anno-db install-feature horse def -name horse_custom
```

For `install-feature`, omitting the version is equivalent to using `def`.
The `-name` mapping is written to `custom.name.tsv` in the selected data
directory.

`-name` also works with `install-genebed` and stores the selected GeneBEDs
under the custom species directory.

`install-feature -o DIR` reuses preprocessed feature files from `DIR`
when `order.lst` and the expected FeatureBED files are already present. Like
`download-genebed`, feature generation reuses the cached GTF, then an expected GTF
from the output directory or current working directory before downloading.
`install-feature` reuses existing generated files by default; pass
`--overwrite` to rebuild the installed annotation BED prerequisite and feature
files.

For UCSC short genome IDs, select the UCSC gene table when needed:

```bash
# Use the UCSC RefSeq gene annotation for the T2T human build.
sjcab-peak2anno-db download-genebed hs1 1 --ucsc-source refseq
# Use the UCSC RefSeq gene annotation for the alpaca build.
sjcab-peak2anno-db download-genebed vicPac2 1 --ucsc-source refseq
# Use the UCSC Ensembl gene annotation for the yeast build.
sjcab-peak2anno-db download-genebed sacCer3 R64-1-1 --ucsc-source ens
```

UCSC files are selected from `bigZips/genes/`; `ens` uses `ensGene` and
`refseq` uses `ncbiRefSeq` or `refGene`. GENCODE is preferred when the genome
has a configured GENCODE URL; otherwise recognized UCSC builds are tried before
Ensembl Vertebrates and Ensembl Genomes. See the
[UCSC downloads page](https://hgdownload.soe.ucsc.edu/downloads.html).

Default installed layout:

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
- `--promoter-down`: downstream promoter flank; defaults to `--promoter-bp`.
- `--distal-bp`/`-D`: distal flank size, default `50kb`.
- `--distal-down`: downstream distal flank; defaults to `--distal-bp`.
- `--tes-bp`/`-e`: TES flank size, default `2kb`.
- `--tes-up`: upstream TES flank; defaults to `--tes-bp`.
- `--prefix`/`-P`: output prefix, default is the promoter size label.
- `--gene-bed`/`-b`: use an existing GeneBED for region generation.
- `--gtf-path`/`-g`: use an existing local GTF.

Python API:

```python
import sjcab_peak2anno_db as db

db.download_feature("hg38", "v31", "feature_downloads")

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

```bash
# Install Roadmap ChromHMM resources for selected genomes and samples.
sjcab-peak2anno-db install-chromhmm -m 18 -s hg19 -i E001,E063
sjcab-peak2anno-db install-chromhmm -m 18 -s hg38 -i E063

sjcab-peak2anno-db download-chromhmm -o chromhmm_downloads -m 25 -s hg19 -c GM12878
sjcab-peak2anno-db download-chromhmm -o chromhmm_downloads -m 18 -s hg38 -i E063
sjcab-peak2anno-db download-chromhmm -o chromhmm_downloads -m 15 -t brain
```

ChromHMM defaults to Roadmap hg19 BED files. Use `--species hg38`/`-s hg38` for
Roadmap lifted-over BEDs. Other UCSC builds use `--yes-liftover` and the hg38
chain files listed at `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/`.
Select records with epigenome IDs (`--ids`/`-i`),
fuzzy tissue/group text (`--tissue`/`-t`), or fuzzy sample/cell-line text
(`--cellline`/`-c`).

Output layout:

```text
{db-path}/chromhmm/{genome}/{model}state/{epigenome_id}_{model}_{model_label}_dense.bed.gz
{db-path}/chromhmm/metadata.tsv
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

With no selector, the command downloads
`segway_encyclopedia.bed.gz`. Existing files are skipped by default; pass
`--overwrite` to replace them.

```bash
# Install or download Segway annotation resources.
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
sjcab-peak2anno-db download-segway -o segway_downloads -s hg38 --yes-liftover -i GM12878
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
{db-path}/segway/hg19/segway_encyclopedia.bed.gz
{db-path}/segway/hg19/caas.bed.gz
{db-path}/segway/hg19/label_info.tab
{db-path}/segway/hg19/{sample}.bed.gz
{db-path}/segway/metadata.tsv
{db-path}/segway/liftover_hg19_to_{genome}.sh
```

Lifted BED files include `lift` in the genome suffix, for example:

```text
{db-path}/segway/hg38/PERIPHERAL_BLOOD_MONONUCLEAR_PRIMARY_CELLS.hg38lift.bed.gz
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
{db-path}/blacklists/{name}.bed.20230411
{db-path}/blacklists/{name}.bed
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
{db-path}/cgi/{species}_cgi.bed
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
sjcab-peak2anno-db install genebed feature blacklists cgi
sjcab-peak2anno-db install -c genebed -c feature
sjcab-peak2anno-db install --overwrite

sjcab-peak2anno-db install-genebed
sjcab-peak2anno-db download-genebed hg38 v31 -o annotations
sjcab-peak2anno-db download-genebed mm10 vM22 -g gencode.vM22.annotation.gtf.gz -o annotations
sjcab-peak2anno-db dedup-bed hg38 v31 -m peak -i h3k4me3_peaks.bed -o annotations
sjcab-peak2anno-db dedup-bed mm10 vM22 -m isoID -i isoforms.txt -K ensid -o annotations
sjcab-peak2anno-db filter-bed -b annotations/bed/hg38/v31/all.gene.bed -m perover -i active_chromhmm.bed -o annotations
sjcab-peak2anno-db filter-bed hg38 v31 -m isoexp -i isoform_expression.tsv --exclusive -o annotations

sjcab-peak2anno-db install-feature all all
sjcab-peak2anno-db install-feature hg38 v31 -o feature_downloads
sjcab-peak2anno-db download-feature hg38 v31 -o feature_downloads
sjcab-peak2anno-db download-feature all all -o feature_downloads

sjcab-peak2anno-db install-chromhmm -m 18 -s hg19 -i E001,E063
sjcab-peak2anno-db install-chromhmm -m 18 -s hg38 -i E063
sjcab-peak2anno-db download-chromhmm -o chromhmm_downloads -m 25 -s hg19 -c GM12878
sjcab-peak2anno-db download-chromhmm -o chromhmm_downloads -m 18 -s hg38 -i E063

sjcab-peak2anno-db install-segway -i GM12878,H1-HESC
sjcab-peak2anno-db download-segway -o segway_downloads -t brain
sjcab-peak2anno-db download-segway -o segway_downloads -s hg38 --yes-liftover -i GM12878

sjcab-peak2anno-db install-blacklists
sjcab-peak2anno-db install-cgi
sjcab-peak2anno-db download-cgi
sjcab-peak2anno-db download-cgi -d cgi_downloads --species hg38

sjcab-peak2anno-db path hg38 gene
```
