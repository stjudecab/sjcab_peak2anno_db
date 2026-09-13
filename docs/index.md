# sjcab_peak2anno_db

Gene-backed peak-to-annotation BED resources for SJ/CAB workflows. The package
provides bundled annotations and downloads or generates resources from GENCODE,
UCSC, Ensembl, Roadmap ChromHMM, Segway, blacklist, and CpG-island sources.

Use the navigation to find installation, configuration, resource-generation,
and transcript-selection guides.

## The install/download convention

`install-*` commands put resources under the configured database directory.
`download-*` commands write to a staging directory, which defaults to the
current working directory. Commands ending in `-feature` generate both gene
GeneBED and derived FeatureBED resources.

## Quick start

```bash
pip install sjcab_peak2anno_db
sjcab-peak2anno-db install
sjcab-peak2anno-db list
```

The default database directory is `~/.sjcab_peak2anno_db`. Use `--db-path` or
`SJCAB_PEAK2ANNO_DB_PATH` to select another location.

The original long-form references remain available in the repository:
[README.md](https://github.com/stjudecab/sjcab_peak2anno_db/blob/main/README.md)
and [README.DEDUP.md](https://github.com/stjudecab/sjcab_peak2anno_db/blob/main/README.DEDUP.md).
