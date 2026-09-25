# Segway

Segway downloads the hg19 encyclopedia and can create a liftOver script for
another UCSC assembly. Generated scripts try the UCSC `liftOver` binary from
the `bioconda::ucsc-liftover` package first and fall back to CrossMap.

```bash
sjcab-peak2anno-db install-segway
sjcab-peak2anno-db download-segway -o segway_downloads -i GM12878
sjcab-peak2anno-db download-segway -o segway_downloads -s hg38 --yes-liftover
```

Lifted files and the helper script are stored under the selected database
directory. The script reuses an existing `segway-liftover` environment when
available and caches chain files under `cache/chains`.

## Storage layout

With the default database path, Segway files are stored as:

```text
{db-path}/segway/{genome}/segway_encyclopedia.bed.gz
{db-path}/segway/{genome}/caas.bed.gz
{db-path}/segway/{genome}/{sample}.bed.gz
{db-path}/segway/{genome}/metadata.tsv
{db-path}/segway/segway_downloads.tsv
```

The default `{db-path}` is `~/.sjcab_peak2anno_db`; use `--db-path` or
`SJCAB_PEAK2ANNO_DB_PATH` to change it. LiftOver output is written under
`{db-path}/segway/{target-genome}/`, and generated scripts are stored as
`{db-path}/segway/liftover_hg19_to_{target-genome}.sh`.

## Resource URLs

- [ENCODE Segway software and download page](https://www.encodeproject.org/software/segway/)
- [Segway Encyclopedia publication](https://www.encodeproject.org/publications/94941f71-80c8-43d2-809b-25161efc3be0/)
- [Segway fallback download site](https://noble.gs.washington.edu/proj/encyclopedia/)
- [UCSC hg19 liftOver chains](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/)
