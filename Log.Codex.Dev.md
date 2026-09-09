# Codex Development Log

## 2026-09-08

- Updated `sjcab-peak2anno-db install` with `--ucsc-source` and `--clean-cache`.
- Threaded cache directory, UCSC source, and cache cleanup through the install pipeline.
- Reinstall the package in the target environment before testing the installed console script:
  `/research/rgs01/applications/hpcf/authorized_apps/cab/sjcab_share/condaenv/test/bin/python -m pip install -e .`
