# Command reference

This page lists command options and their defaults. `install-*` commands write
to `--db-path`; `download-*` commands write to `--output-dir`.

## Common options

<div style="overflow-x:auto">
<table>
<colgroup><col style="min-width: 330px"><col style="min-width: 170px"><col style="min-width: 480px"></colgroup>
<thead><tr><th>Option</th><th>Default</th><th>Description</th></tr></thead>
<tbody>
<tr><td style="white-space:nowrap"><code>SPECIES</code>, <code>-s</code>, <code>--species</code></td><td>command-specific</td><td>Species or genome assembly.</td></tr>
<tr><td style="white-space:nowrap"><code>VERSION</code>, <code>-v</code>, <code>--ver</code></td><td><code>def</code> where supported</td><td>Release or version.</td></tr>
<tr><td style="white-space:nowrap"><code>-d</code>, <code>--db-path DIR</code></td><td><code>~/.sjcab_peak2anno_db</code></td><td>Database/cache directory.</td></tr>
<tr><td style="white-space:nowrap"><code>-o</code>, <code>--output-dir DIR</code></td><td><code>.</code> for download commands</td><td>Download or generated-file directory.</td></tr>
<tr><td style="white-space:nowrap"><code>-name NAME</code></td><td>unset</td><td>Custom installed species name.</td></tr>
<tr><td style="white-space:nowrap"><code>-n</code>, <code>--workers N</code></td><td><code>2</code> for selector commands</td><td>Parallel independent TSS/TES output writers.</td></tr>
<tr><td style="white-space:nowrap"><code>-j</code>, <code>--processes N</code></td><td><code>4</code></td><td>Independent GTF workers for install/download feature commands.</td></tr>
<tr><td style="white-space:nowrap"><code>-g</code>, <code>--gtf-path FILE</code></td><td>unset</td><td>Use an existing local GTF.</td></tr>
<tr><td style="white-space:nowrap"><code>-u</code>, <code>--url URL</code></td><td>resolved URL</td><td>Override the GTF URL.</td></tr>
<tr><td style="white-space:nowrap"><code>--ucsc-source {ens,refseq}</code></td><td><code>ens</code></td><td>UCSC gene-table source.</td></tr>
<tr><td style="white-space:nowrap"><code>--clean-cache [DAYS]</code></td><td><code>90</code> when configured</td><td>Remove old cache files after generation.</td></tr>
<tr><td style="white-space:nowrap"><code>--sizes-clean [0|1]</code></td><td><code>1</code></td><td>Create <code>.sizes.clean</code> files.</td></tr>
<tr><td style="white-space:nowrap"><code>-dry-run</code></td><td>off</td><td>Resolve URLs without downloading where supported.</td></tr>
</tbody>
</table>
</div>

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

## `dedup-genebed`, `filter-genebed`, `dedup-feature`, and `filter-feature`

| Option | Default |
| --- | --- |
| `--method`, `-m` | `longcol5` for `dedup-genebed`; required for `filter-genebed` |
| `--promoter-bp`, `-p` | `2kb` |
| `--promoter-down` | same as `--promoter-bp` |
| `--gene-key`, `-K` | `symbol` |
| `--workers`, `-n` | `2` |
| `--exclusive` | off; inclusive matching is used |
| `--output-dir`, `-o` | `.` |

`dedup-feature` and `filter-feature` accept the same selector options as their
`*-bed` counterparts. They retain transcript IDs using the GeneBED selection
and regenerate FeatureBEDs from only those transcripts. `dedup-feature` falls
back to the longest transcript; `filter-feature` omits genes without selector
support.

## Other commands

- `list def`: lists the default bundled GeneBED resources. Use `list genebed`,
  `list feature`, `list blacklists`, `list cgi`, `list chromhmm`, or
  `list segway` to select another group; output is aligned like `column -t`
  and uses a `source` column.
- `path`: prints a resource path; `--install`/`-I` installs missing resources.
- `install-genebed`/`download-genebed`: generate GeneBEDs from GTF input.
- `install-blacklists`, `install-cgi`, `install-chromhmm`, and
  `install-segway`: install the corresponding external resource sets.

Run `sjcab-peak2anno-db COMMAND -h` for command-specific help.
