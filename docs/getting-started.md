# Installation and quick start

## Install the package

```bash
pip install sjcab_peak2anno_db
conda install stjudecab::sjcab_peak2anno_db
```

Install the default resource bundle:

```bash
sjcab-peak2anno-db install
sjcab-peak2anno-db install def
sjcab-peak2anno-db install default
```

Install selected components:

```bash
sjcab-peak2anno-db install genebed feature blacklists cgi
sjcab-peak2anno-db install -c genebed -c feature
```

The supported install components are `genebed`, `feature`, blacklists, and
`cgi`. The three commands above install the complete default bundle. `list` shows
the bundled registry, while `path` prints a resource path.

## Download/install supported species/versions

```bash
sjcab-peak2anno-db download-genebed hg38 v31 -o annotations
sjcab-peak2anno-db download-feature mm10 vM22 -o annotations
sjcab-peak2anno-db install-feature cat def
sjcab-peak2anno-db install-feature dog 115
sjcab-peak2anno-db install-feature vicpac2
```

Use `-j/--processes` to process independent GTFs in parallel. GTF downloads are
sequential with a short randomized pause, and conversion within each worker is
single-process. The worker count is capped at the CPUs allocated to the job:

```bash
sjcab-peak2anno-db download-feature cat,dog,horse -j 3
```

`dedup-genebed` and `filter-genebed` support `-n/--workers` (default `2`) to write
their independent TSS and TES outputs in parallel.

## Python API

```python
import sjcab_peak2anno_db as db

db.install_data()
db.download_and_convert_gencode_gtf("hg38", "v31", output_dir="annotations")
db.download_feature("hg38", "v31", "annotations")
```

## Data lookup and cache behavior

Downloads first search the database `{db-path}/cache`, then the output directory, then the
current working directory. Download URLs are recorded in
`{db-path}/download_urls.log`. Match of gtf file name and species/version listed in 
`{db-path}/installed.tsv`. GENCODE, UCSC, and Ensembl resolution details
are documented in the [GeneBEDs](genebeds.md) guide.
