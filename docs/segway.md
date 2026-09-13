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

## Resource URLs

- [ENCODE Segway software and download page](https://www.encodeproject.org/software/segway/)
- [Segway Encyclopedia publication](https://www.encodeproject.org/publications/94941f71-80c8-43d2-809b-25161efc3be0/)
- [Segway fallback download site](https://noble.gs.washington.edu/proj/encyclopedia/)
- [UCSC hg19 liftOver chains](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/)
