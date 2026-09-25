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

The generated database-level helper is `{db-path}/liftover_hg38_to.sh`:

```bash
bash liftover_hg38_to.sh
bash liftover_hg38_to.sh mm10 crossmap
bash liftover_hg38_to.sh mm10 ucsc
```

It scans the blacklist, CGI, ChromHMM, and Segway directories and skips files
that already have lifted outputs. CrossMap is the default; `auto` tries
CrossMap and then UCSC `liftOver`.
The helper detects package managers in this order: Pixi, micromamba, mamba,
then conda, and installs the required tool when it is unavailable.

Output is stored under `{db-path}/chromhmm/{genome}` with a `metadata.tsv`
index.

## Resource URLs

- [Roadmap metadata](https://egg2.wustl.edu/roadmap/web_portal/data.js)
- [Roadmap ChromHMM core 15-state models](https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#core_15state)
- [Roadmap 15-state BED files](https://egg2.wustl.edu/roadmap/data/byFileType/)
- [UCSC hg38 liftOver chains](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/)
