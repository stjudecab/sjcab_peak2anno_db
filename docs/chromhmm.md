# ChromHMM

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

Output is stored under `{db-path}/chromhmm/{genome}` with a `metadata.tsv`
index.
