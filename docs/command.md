# Command reference

This page lists command options and their defaults. `install-*` commands write
to `--db-path`; `download-*` commands write to `--output-dir`.

## Common options

| Option | Default | Description |
| --- | --- | --- |
| `SPECIES`, `-s`, `--species` | command-specific | Species or genome assembly. |
| `VERSION`, `-v`, `--ver` | `def` where supported | Release or version. |
| `-d`, `--db-path DIR` | `~/.sjcab_peak2anno_db` | Database/cache directory. |
| `-o`, `--output-dir DIR` | `.` for download commands | Download or generated-file directory. |
| `-name NAME` | unset | Custom installed species name. |
| `-n`, `--workers N` | `2` for `dedup-bed`/`filter-bed` | Parallel independent TSS/TES output writers. |
| `-j`, `--processes N` | `4` | Independent GTF workers for install/download feature commands. |
| `-g`, `--gtf-path FILE` | unset | Use an existing local GTF. |
| `-u`, `--url URL` | resolved URL | Override the GTF URL. |
| `--ucsc-source {ens,refseq}` | `ens` | UCSC gene-table source. |
| `--clean-cache [DAYS]` | `90` when configured | Remove old cache files after generation. |
| `--sizes-clean [0|1]` | `1` | Create `.sizes.clean` files. |
| `-dry-run` | off | Resolve URLs without downloading where supported. |

`-n/--workers` is intentionally separate from `-j/--processes`: workers apply
to the independent left/right-style TSS/TES outputs produced by the selector
commands, while processes apply to independent GTF jobs.

## `install`

```text
sjcab-peak2anno-db install [COMPONENT ...]
```

Components are `genebed`, `feature`, `blacklists`, and `cgi`. With no component,
or with `def` or `default`, all configured default components are installed.
The default overwrite behavior is off; use `--overwrite` to regenerate files.

## `install-feature` and `download-feature`

| Option | Default |
| --- | --- |
| `--promoter-bp`, `-p` | `2kb` |
| `--promoter-down` | same as `--promoter-bp` |
| `--distal-bp`, `-D` | `50kb` |
| `--distal-down` | same as `--distal-bp` |
| `--tes-bp`, `-e` | `2kb` |
| `--tes-up` | same as `--tes-bp` |
| `--prefix`, `-P` | promoter-size label |
| `--include-tss-base`, `-B` | off |
| `--no-overwrite`, `-n` | off; existing files may be regenerated |
| `--processes`, `-j` | `4` |

## `dedup-bed` and `filter-bed`

| Option | Default |
| --- | --- |
| `--method`, `-m` | `longcol5` for `dedup-bed`; required for `filter-bed` |
| `--promoter-bp`, `-p` | `2kb` |
| `--promoter-down` | same as `--promoter-bp` |
| `--gene-key`, `-K` | `symbol` |
| `--workers`, `-n` | `2` |
| `--exclusive` | off; inclusive matching is used |
| `--output-dir`, `-o` | `.` |

## Other commands

- `list def`: lists the default bundled GeneBED resources. Use `list genebed`,
  `list feature`, `list blacklists`, or `list cgi` to select another group;
  output is aligned like `column -t` and uses a `source` column.
- `path`: prints a resource path; `--install`/`-I` installs missing resources.
- `install-genebed`/`download-genebed`: generate GeneBEDs from GTF input.
- `install-blacklists`, `install-cgi`, `install-chromhmm`, and
  `install-segway`: install the corresponding external resource sets.

Run `sjcab-peak2anno-db COMMAND -h` for command-specific help.
