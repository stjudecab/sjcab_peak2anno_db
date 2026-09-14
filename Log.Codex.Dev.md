# Codex Development Log

## 2026-09-08

- Updated `sjcab-peak2anno-db install` with `--ucsc-source` and `--clean-cache`.
- Threaded cache directory, UCSC source, and cache cleanup through the install pipeline.
- Switched UCSC chromosome-size generation to NCBI assembly reports.
- Added matching between UCSC builds, NCBI assemblies, and organism metadata.
- Added fallbacks for empty clean-size files and refreshed packaged Ensembl
  clean-size data, removing obsolete `fr2` data.
- Added retry handling for Ensembl chromosome metadata requests and taxonomy-
  based Ensembl assembly resolution.
- Added the initial RC/configuration and annotation-resolution improvements
  captured by the checkpoint commits for this day.

## 2026-09-09

- Added Ensembl cache refresh logic and feature URL options.
- Improved GTF, assembly, and chromosome-size resolution behavior.

## 2026-09-10

- Added RC-file configuration with XDG and home-file lookup, automatic creation
  of a commented template, environment-variable overrides, install-component
  defaults, species/version pairing, stale-cache controls, size-clean controls,
  and cache cleanup policies.
- Added liftOver support and chain-file caching for external resource types.
- Removed packaged hg38 rDNA blacklist files and updated blacklist handling.
- Added support for custom GeneBED/FeatureBED naming, assembly-based resource
  selection, and improved CLI configuration handling.

## 2026-09-11

- Reorganized README command options and common configuration guidance.
- Released versions `0.1.7` and `0.1.8`.
- Added and updated package publishing workflow support, including manual
  workflow dispatch.

## 2026-09-12

- Added parallel GTF-to-GeneBED and FeatureBED generation with configurable
  worker processes.
- Added MkDocs and Read the Docs configuration and reorganized documentation
  into separate topic pages.
- Standardized terminology from “genebed”/“feature bed” to GeneBEDs and
  FeatureBEDs throughout the CLI and documentation.
- Added directional promoter, TES, and distal-window controls and custom naming
  for GeneBEDs and FeatureBEDs.
- Added configurable BED/text score columns and peak-text delimiter handling.
- Added resource documentation and links for blacklists, CGI, ChromHMM,
  Segway, UCSC, ENCODE, Roadmap, and liftOver chains.
- Removed explicit delimiter lists from text-selector guidance and documented
  environment-variable defaults.

## 2026-09-13

- Changed `peak` and `perover` selectors to detect BED versus text once per
  file from the first data row, then apply that mode consistently to all rows.
- Updated selector documentation to describe per-file detection.

## 2026-09-14

- Added automatic VCS-based development version metadata.
