# ChromHMM, Segway, and external resources

## ChromHMM

ChromHMM commands download Roadmap dense-state BED files and normalized
metadata. Directly supported species include `hg19`, `hg38`, `mm39`, `mm10`,
and `mm9`; other UCSC builds can use liftOver from hg38.

```bash
sjcab-peak2anno-db install-chromhmm -m 18 -s hg19 -i E001,E063
sjcab-peak2anno-db download-chromhmm -o chromhmm_downloads -m 25 -s hg19 -c GM12878
sjcab-peak2anno-db download-chromhmm -o chromhmm_downloads -m 18 -s hg38 -i E063
```

Use `--ids`/`-i`, `--tissue`/`-t`, and `--cellline`/`-c` to select records.
Other UCSC builds use `--yes-liftover`; chain files are cached under
`{db-path}/cache/chains/`.

Output is stored under `{db-path}/chromhmm/{genome}` with a
`metadata.tsv` index.

## Segway

Segway downloads the hg19 encyclopedia and can create a CrossMap liftOver
script for another UCSC assembly.

```bash
sjcab-peak2anno-db install-segway
sjcab-peak2anno-db download-segway -o segway_downloads -i GM12878
sjcab-peak2anno-db download-segway -o segway_downloads -s hg38 --yes-liftover
```

Lifted files and the helper script are stored under the selected database
directory. The script reuses an existing `segway-liftover` environment when
available and caches chain files under `cache/chains`.

## Blacklists and CGI

Install bundled blacklists and CpG-island BEDs:

```bash
sjcab-peak2anno-db install-blacklists
sjcab-peak2anno-db install-cgi
```

Download fresh UCSC CGI tables:

```bash
sjcab-peak2anno-db download-cgi --species hg38 -d cgi_downloads
```

Outputs are stored as `{db-path}/blacklists/{name}.bed` and
`{db-path}/cgi/{species}_cgi.bed`. For blacklist provenance, see the
[bundled blacklist notes](https://github.com/stjudecab/sjcab_peak2anno_db/blob/main/src/sjcab_peak2anno_db/data/blacklists/readme.md).
