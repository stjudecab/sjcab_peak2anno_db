# Common options and configuration

## Common command options

- `SPECIES`, `-s`, `--species`: species name or genome assembly. Matching uses
  GENCODE, then UCSC, then Ensembl. Ensembl candidates are printed when a
  unique match is not found.
- `VERSION`, `-v`, `--ver`: release or version, including `def`, `default`,
  `current`, and `latest` where supported.
- `-d`, `--db-path`: database and cache directory.
- `-o`, `--output-dir`: staging directory for download commands.
- `-name NAME`: store selected GeneBEDs or FeatureBEDs under a custom name;
  records the mapping in `custom.name.tsv`.
- `-n`, `--no-overwrite`: preserve existing generated files.
- `--clean-cache [DAYS]`: remove old cache files after generation. The default
  is 90 days; a negative value removes the current cached GTF immediately.
- `--sizes-clean [0|1]`: create `.sizes.clean` files by default; use `0` to
  disable them.
- `-j`, `--processes N`: worker processes, with one GTF handled by each worker.
  GTF downloads are sequential with a short randomized pause, and conversion
  within each worker is single-process. The default is `4`, capped at the CPUs
  allocated to the current job.
- `-n`, `--workers N`: workers for independent TSS/TES outputs from
  `dedup-bed` and `filter-bed`. The default is `2`.
- `-g`, `--gtf-path`: use a local GTF.
- `-u`, `--url`: override the resolved GTF URL.
- `--ucsc-source {ens,refseq}`: choose the UCSC gene table.
- `-dry-run`: resolve and print URLs without downloading where supported.

Species and versions accept comma-separated values or a `.lst`/`.list` file.
Two-column files and `species:version` values are supported. A single value is
broadcast across the other list; equal-length lists are paired in order.

## RC files and environment variables

The generated configuration template is searched in this order:

1. `SJCAB_PEAK2ANNO_CONFIG`, when set;
2. `$XDG_CONFIG_HOME/sjcab_peak2anno/.sjcab_peak2anno.rc`;
3. `~/.sjcab_peak2anno.rc`.

If the XDG file does not exist, a fully commented template is created. Existing
RC files are completed with any missing commented settings. Commented
configuration variables use `#SJCAB_...` with no space after `#`. RC keys use
the same names as environment variables. The environment prefix is
`SJCAB_PEAK2ANNO_DB_`.

### Environment variable defaults are below can be override

```text
SJCAB_PEAK2ANNO_DB_PATH=~/.sjcab_peak2anno_db
SJCAB_PEAK2ANNO_DB_INSTALL_OPTIONS=genebed,feature,blacklists,cgi
SJCAB_PEAK2ANNO_DB_INSTALL_SPECIES=hg38,hg19,mm10,mm39
SJCAB_PEAK2ANNO_DB_INSTALL_VERSIONS=v31,v31lift37,vM22,vM39
SJCAB_PEAK2ANNO_DB_VERSION_STALE_DAYS=90
SJCAB_PEAK2ANNO_DB_SIZESCLEAN=1
SJCAB_PEAK2ANNO_DB_CLEANCACHE=90
SJCAB_PEAK2ANNO_DB_PEAK_TXT_DELIMITER=:-*/^;_%$,
SJCAB_PEAK2ANNO_DB_BED_SCORE_COLUMN=5
SJCAB_PEAK2ANNO_DB_TXT_SCORE_COLUMN=2
```

`VERSION_STALE_DAYS` controls refresh of Ensembl `VERSION` and `SPECIES`
catalogues. `SIZESCLEAN=0` disables `.sizes.clean`; `CLEANCACHE` controls cache
cleanup age. Explicit command-line options take precedence.

`PEAK_TXT_DELIMITER` controls delimiters inside text-mode `peak` regions; it is
not a command-line option. Selector format is auto-detected once per file from
the first data none-header row: BED selectors use the configured `BED_SCORE_COLUMN`
(default `5`), while text selectors use `TXT_SCORE_COLUMN` (default `2`).
