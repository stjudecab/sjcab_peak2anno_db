# Segway

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
