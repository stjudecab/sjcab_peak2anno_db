# CGI

Install packaged CpG-island BEDs or download fresh UCSC `cpgIslandExt`
tables:

```bash
sjcab-peak2anno-db install-cgi
sjcab-peak2anno-db install cgi
sjcab-peak2anno-db download-cgi --species hg38 -d cgi_downloads
```

Outputs are stored as `{db-path}/cgi/{species}_cgi.bed`.

## Resource URLs

- [UCSC genome downloads](https://hgdownload.soe.ucsc.edu/downloads.html)
- [UCSC per-assembly database tables](https://hgdownload.soe.ucsc.edu/goldenPath/)
- [UCSC CpG island table pattern](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz)
