# Installation and quick start

## Install the package

```bash
pip install sjcab_peak2anno_db
conda install stjudecab::sjcab_peak2anno_db
```

Install the default resource bundle:

```bash
sjcab-peak2anno-db install
```

Install selected components:

```bash
sjcab-peak2anno-db install genebed feature blacklists cgi
sjcab-peak2anno-db install -c genebed -c feature
```

The supported install components are `genebed`, `feature`, `blacklists`, and
`cgi`. `list` shows the bundled registry, while `path` prints a resource path.

## First downloaded resources

```bash
sjcab-peak2anno-db download-genebed hg38 v31 -o annotations
sjcab-peak2anno-db download-feature mm10 vM22 -o annotations
sjcab-peak2anno-db install-feature cat def
```

Use `-j/--processes` to parallelize GTF-to-genebed conversion and feature BED
generation:

```bash
sjcab-peak2anno-db download-feature hg38 v31 -j 8
```

## Python API

```python
import sjcab_peak2anno_db as db

db.install_data()
db.download_and_convert_gencode_gtf("hg38", "v31", output_dir="annotations")
db.download_feature("hg38", "v31", "annotations")
```

## Data lookup and cache behavior

Downloads first search the database cache, then the output directory, then the
current working directory. Download URLs are recorded in
`{db-path}/download_urls.log`. GENCODE, UCSC, and Ensembl resolution details
are documented in the [Genebed](genebed.md) guide.
