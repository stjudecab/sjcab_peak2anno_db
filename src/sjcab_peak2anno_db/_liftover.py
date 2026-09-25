"""Generate UCSC liftOver/CrossMap helpers for BED resources."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional, Union

PathLike = Union[str, Path]
UCSC_LIFTOVER_ROOT = "https://hgdownload.soe.ucsc.edu/goldenPath/{}/liftOver/"


def liftover_script(
    input_dir: PathLike,
    output_dir: PathLike,
    source_genome: str,
    target_genome: str,
    env_name: str = ".sjcab-liftover",
    install_root: Optional[PathLike] = None,
    rename_genome_prefix: bool = False,
    cache_dir: Optional[PathLike] = None,
    mode: str = "crossmap",
) -> str:
    """Return a bash script that lifts BED and BED.GZ files in a directory."""

    source = _genome(source_genome)
    target = _genome(target_genome)
    if source == target:
        raise ValueError("Source and target genomes must differ for liftover.")
    selected_mode = str(mode).strip().lower()
    if selected_mode not in {"crossmap", "ucsc", "auto"}:
        raise ValueError("mode must be crossmap, ucsc, or auto")
    chain_name = "{}To{}.over.chain.gz".format(source, _ucsc_token(target))
    from ._gencode import ucsc_liftover_chain_url
    from ._registry import user_data_dir

    chain_url = ucsc_liftover_chain_url(source, chain_name, cache_dir=cache_dir)
    chain_cache_dir = user_data_dir(cache_dir) / "cache" / "chains"
    install_root_text = (
        _shell_quote(str(Path(install_root).expanduser()))
        if install_root is not None
        else "''"
    )
    output_root = Path(output_dir).expanduser().parent
    return """#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob globstar

ENV_NAME={env_name}
DEFAULT_TARGET={target}
TARGET="${{1:-${{DEFAULT_TARGET}}}}"
REQUESTED_MODE="${{2:-{mode}}}"
INPUT_DIR={input_dir}
OUTPUT_ROOT={output_root}
INSTALL_ROOT={install_root}
CHAIN_DIR={chain_dir}
TMP_DIR="${{TMPDIR:-/lustre_scratch/user_scratch/bxu2/TMPDIR/codex}}/sjcab-liftover.$$"
mkdir -p "${{TMP_DIR}}"
trap 'rm -rf "${{TMP_DIR}}"' EXIT

case "${{REQUESTED_MODE,,}}" in
  crossmap|ucsc|auto) ;;
  *) echo "Mode must be crossmap, ucsc, or auto: ${{REQUESTED_MODE}}" >&2; exit 2 ;;
esac

ucsc_token() {{
  local value="$1" lower
  lower="${{value,,}}"
  case "${{lower}}" in
    canfam*) echo "CanFam${{value:6}}" ;;
    equcab*) echo "EquCab${{value:6}}" ;;
    bostau*) echo "BosTau${{value:6}}" ;;
    galgal*) echo "GalGal${{value:6}}" ;;
    xenla*) echo "XenLa${{value:5}}" ;;
    danrer*) echo "DanRer${{value:6}}" ;;
    susscr*) echo "SusScr${{value:6}}" ;;
    saccer*) echo "SacCer${{value:6}}" ;;
    *) echo "${{value^}}" ;;
  esac
}}

CHAIN_NAME="{source}To$(ucsc_token "${{TARGET}}").over.chain.gz"
CHAIN="${{CHAIN_DIR}}/${{CHAIN_NAME}}"
CHAIN_URL="{chain_root}${{CHAIN_NAME}}"
OUTPUT_DIR="${{OUTPUT_ROOT}}/${{TARGET}}"
mkdir -p "${{CHAIN_DIR}}" "${{OUTPUT_DIR}}"
if [ ! -s "${{CHAIN}}" ]; then
  curl -L "${{CHAIN_URL}}" -o "${{CHAIN}}"
fi

activate_existing_env() {{
  if command -v micromamba >/dev/null 2>&1; then
    hook="$(micromamba shell hook -s bash 2>/dev/null)" && eval "${{hook}}"
    micromamba activate "${{ENV_NAME}}" >/dev/null 2>&1 && return 0
  fi
  if command -v mamba >/dev/null 2>&1; then
    hook="$(mamba shell hook -s bash 2>/dev/null)" && eval "${{hook}}"
    mamba activate "${{ENV_NAME}}" >/dev/null 2>&1 && return 0
  fi
  if command -v conda >/dev/null 2>&1; then
    hook="$(conda shell.bash hook 2>/dev/null)" && eval "${{hook}}"
    conda activate "${{ENV_NAME}}" >/dev/null 2>&1 && return 0
  fi
  return 1
}}

create_and_activate_env() {{
  if command -v micromamba >/dev/null 2>&1; then
    if micromamba create -y -n "${{ENV_NAME}}" -c conda-forge -c bioconda crossmap; then
      eval "$(micromamba shell hook -s bash)"
      micromamba activate "${{ENV_NAME}}"
      return 0
    fi
    micromamba create -y -n "${{ENV_NAME}}" bioconda::ucsc-liftover
    eval "$(micromamba shell hook -s bash)"
    micromamba activate "${{ENV_NAME}}"
  elif command -v mamba >/dev/null 2>&1; then
    if mamba create -y -n "${{ENV_NAME}}" -c conda-forge -c bioconda crossmap; then
      eval "$(mamba shell hook -s bash)"
      mamba activate "${{ENV_NAME}}"
      return 0
    fi
    mamba create -y -n "${{ENV_NAME}}" bioconda::ucsc-liftover
    eval "$(mamba shell hook -s bash)"
    mamba activate "${{ENV_NAME}}"
  elif command -v conda >/dev/null 2>&1; then
    if conda create -y -n "${{ENV_NAME}}" -c conda-forge -c bioconda crossmap; then
      eval "$(conda shell.bash hook)"
      conda activate "${{ENV_NAME}}"
      return 0
    fi
    conda create -y -n "${{ENV_NAME}}" bioconda::ucsc-liftover
    eval "$(conda shell.bash hook)"
    conda activate "${{ENV_NAME}}"
  else
    echo "micromamba, mamba, or conda is required to install liftOver or CrossMap" >&2
    exit 1
  fi
}}

if command -v pixi >/dev/null 2>&1; then
  PIXI_MODE=1
elif ! activate_existing_env; then
  create_and_activate_env
fi

if [ "${{PIXI_MODE:-0}}" = "1" ]; then
  LIFTOVER_MODE=pixi
elif [ "${{REQUESTED_MODE,,}}" = "crossmap" ] || [ "${{REQUESTED_MODE,,}}" = "auto" ]; then
  if command -v CrossMap >/dev/null 2>&1; then
    CROSSMAP=(CrossMap)
    LIFTOVER_MODE=crossmap
  elif command -v CrossMap.py >/dev/null 2>&1; then
    CROSSMAP=(CrossMap.py)
    LIFTOVER_MODE=crossmap
  fi
fi
if [ -z "${{LIFTOVER_MODE:-}}" ] && {{ [ "${{REQUESTED_MODE,,}}" = "ucsc" ] || [ "${{REQUESTED_MODE,,}}" = "auto" ]; }} && command -v liftOver >/dev/null 2>&1; then
  LIFTOVER=(liftOver)
  LIFTOVER_MODE=ucsc
fi
if [ -z "${{LIFTOVER_MODE:-}}" ] && command -v CrossMap >/dev/null 2>&1; then
  CROSSMAP=(CrossMap)
  LIFTOVER_MODE=crossmap
fi
if [ -z "${{LIFTOVER_MODE:-}}" ] && command -v CrossMap.py >/dev/null 2>&1; then
  CROSSMAP=(CrossMap.py)
  LIFTOVER_MODE=crossmap
fi
if [ -z "${{LIFTOVER_MODE:-}}" ]; then
  echo "liftOver and CrossMap are not available in conda environment ${{ENV_NAME}}." >&2
  exit 1
fi

run_crossmap() {{
  local chain="$1" input="$2" output="$3" fields padded mapped
  fields="$(awk 'NF {{ print NF; exit }}' "${{input}}")"
  if [ "${{fields:-0}}" -lt 12 ]; then
    padded="${{TMP_DIR}}/$(basename "${{input}}").bed12"
    mapped="${{TMP_DIR}}/$(basename "${{output}}").mapped"
    awk 'BEGIN {{ OFS="\t" }} NF {{ print $0, ".", 0, "+", $2, $3, 0, 1, $3-$2, 0 }}' "${{input}}" > "${{padded}}"
    if [ "${{LIFTOVER_MODE}}" = "pixi" ]; then
      pixi exec --channel conda-forge --channel bioconda --spec crossmap CrossMap bed "${{chain}}" "${{padded}}" "${{mapped}}"
    else
      "${{CROSSMAP[@]}}" bed "${{chain}}" "${{padded}}" "${{mapped}}"
    fi
    awk -v n="${{fields}}" 'BEGIN {{ OFS="\t" }} NF {{ for (i=1; i<=n; i++) {{ printf "%s%s", $i, (i==n ? ORS : OFS) }} }}' "${{mapped}}" > "${{output}}"
  elif [ "${{LIFTOVER_MODE}}" = "pixi" ]; then
    pixi exec --channel conda-forge --channel bioconda --spec crossmap CrossMap bed "${{chain}}" "${{input}}" "${{output}}"
  else
    "${{CROSSMAP[@]}}" bed "${{chain}}" "${{input}}" "${{output}}"
  fi
}}

run_ucsc() {{
  local chain="$1" input="$2" output="$3" fields padded mapped unmapped original
  fields="$(awk 'NF {{ print NF; exit }}' "${{input}}")"
  if [ "${{fields:-0}}" -lt 12 ]; then
    padded="${{TMP_DIR}}/$(basename "${{input}}").bed12"
    mapped="${{TMP_DIR}}/$(basename "${{output}}").mapped"
    unmapped="${{TMP_DIR}}/$(basename "${{output}}").unmap"
    original="${{TMP_DIR}}/$(basename "${{input}}").original"
    awk -v padded="${{padded}}" -v original="${{original}}" 'BEGIN {{ OFS="\t" }} NF {{
      id="sjcab_liftover_" NR
      print $1, $2, $3, id, 0, "+", $2, $3, 0, 1, $3-$2, 0 > padded
      printf "%s", id > original
      for (i=4; i<=NF; i++) printf "%s%s", OFS, $i > original
      printf "%s", ORS > original
    }}' "${{input}}"
    if [ "${{LIFTOVER_MODE:-}}" = "pixi" ]; then
      pixi exec --channel conda-forge --channel bioconda --spec ucsc-liftover liftOver "${{padded}}" "${{chain}}" "${{mapped}}" "${{unmapped}}"
    else
      "${{LIFTOVER[@]}}" "${{padded}}" "${{chain}}" "${{mapped}}" "${{unmapped}}"
    fi
    for source in "${{mapped}}" "${{unmapped}}"; do
      target="${{output}}"
      [ "${{source}}" = "${{unmapped}}" ] && target="${{output}}.unmap"
      awk -v n="${{fields}}" -v original="${{original}}" 'BEGIN {{
        OFS="\t"
        while ((getline line < original) > 0) {{
          m=split(line, a, "\t")
          for (i=2; i<=m; i++) saved[a[1], i-1]=a[i]
        }}
        close(original)
      }} /^#/ {{ print; next }} NF {{
        printf "%s\t%s\t%s", $1, $2, $3
        for (i=4; i<=n; i++) printf "%s%s", OFS, saved[$4, i]
        print ""
      }}' "${{source}}" > "${{target}}"
    done
  else
    "${{LIFTOVER[@]}}" "${{input}}" "${{chain}}" "${{output}}" "${{output}}.unmap"
  fi
}}

run_liftover() {{
  if [ "${{LIFTOVER_MODE}}" = "pixi" ]; then
    if [ "${{REQUESTED_MODE,,}}" = "ucsc" ]; then
      run_ucsc "$@"
    else
      run_crossmap "$@"
    fi
  elif [ "${{LIFTOVER_MODE}}" = "ucsc" ]; then
    run_ucsc "$@"
  else
    run_crossmap "$@"
  fi
}}

while IFS= read -r -d '' bed; do
  rel="${{bed#${{INPUT_DIR}}/}}"
  if [ "{rename_genome_prefix}" = "1" ]; then
    rel="${{TARGET}}${{rel#{source}}}"
  fi
  out="${{OUTPUT_DIR}}/${{rel}}"
  if [[ "${{out}}" == *.gz ]]; then
    plain="${{out%.gz}}"
    [ -s "${{out}}" ] && continue
  else
    plain="${{out}}"
    [ -s "${{out}}" ] && continue
  fi
  mkdir -p "$(dirname "${{plain}}")"
  if [[ "${{bed}}" == *.gz ]]; then
    unpacked="${{TMP_DIR}}/$(basename "${{plain}}")"
    gzip -dc "${{bed}}" > "${{unpacked}}"
    run_liftover "${{CHAIN}}" "${{unpacked}}" "${{plain}}"
    gzip -f "${{plain}}"
  else
    run_liftover "${{CHAIN}}" "${{bed}}" "${{plain}}"
  fi
  [ -e "${{plain}}.unmap" ] || : > "${{plain}}.unmap"
done < <(find "${{INPUT_DIR}}" -path "${{INPUT_DIR}}/.liftover" -prune -o \\( -type f -o -type l \\) {source_filter} \\( -name '*.bed' -o -name '*.bed.gz' \\) -print0)

if [ -n "${{INSTALL_ROOT}}" ]; then
  mkdir -p "${{INSTALL_ROOT}}"
  cp -a "${{OUTPUT_DIR}}"/. "${{INSTALL_ROOT}}/"
  rm -rf "${{OUTPUT_DIR}}"
  rmdir "$(dirname "${{OUTPUT_DIR}}")" 2>/dev/null || true
fi
""".format(
        chain_name=chain_name,
        chain_root=UCSC_LIFTOVER_ROOT.format(source),
        chain_dir=_shell_quote(str(chain_cache_dir)),
        env_name=_shell_quote(env_name),
        input_dir=_shell_quote(str(Path(input_dir).expanduser())),
        output_root=_shell_quote(str(output_root)),
        install_root=install_root_text,
        output_dir=_shell_quote(str(Path(output_dir).expanduser())),
        source_filter=(
            "-iname {}".format(_shell_quote(source + "*"))
            if rename_genome_prefix
            else ""
        ),
        rename_genome_prefix="1" if rename_genome_prefix else "0",
        source=source,
        target=target,
        mode=selected_mode,
    )


def write_liftover_script(
    input_dir: PathLike,
    output_dir: PathLike,
    source_genome: str,
    target_genome: str,
    script_path: PathLike,
    env_name: str = ".sjcab-liftover",
    install_root: Optional[PathLike] = None,
    rename_genome_prefix: bool = False,
    cache_dir: Optional[PathLike] = None,
    mode: str = "crossmap",
) -> Path:
    """Write :func:`liftover_script` atomically and make it executable."""

    output = Path(script_path).expanduser()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    try:
        temporary.write_text(
            liftover_script(
                input_dir,
                output_dir,
                source_genome,
                target_genome,
                env_name=env_name,
                install_root=install_root,
                rename_genome_prefix=rename_genome_prefix,
                cache_dir=cache_dir,
                mode=mode,
            ),
            encoding="utf-8",
        )
        temporary.chmod(0o755)
        temporary.replace(output)
    finally:
        if temporary.exists():
            temporary.unlink()
    return output


def resource_liftover_script(
    data_dir: PathLike,
    target_genome: str = "hg38",
    script_path: Optional[PathLike] = None,
    mode: str = "crossmap",
) -> Path:
    """Write one helper that scans all external-resource directories.

    Blacklists, CGI, and ChromHMM source files are hg38-based. Segway files
    are hg19-based, so the generated helper selects the matching chain for
    each source tree while writing all outputs into their normal resource
    folders.
    """

    selected_mode = str(mode).strip().lower()
    if selected_mode not in {"crossmap", "ucsc", "auto"}:
        raise ValueError("mode must be crossmap, ucsc, or auto")
    root = Path(data_dir).expanduser()
    output = Path(script_path).expanduser() if script_path else root / "liftover_hg38_to.sh"
    output.parent.mkdir(parents=True, exist_ok=True)
    script = r'''#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob globstar

ROOT={root}
DEFAULT_TARGET={target}
TARGET="${{1:-${{DEFAULT_TARGET}}}}"
REQUESTED_MODE="${{2:-{mode}}}"
CHAIN_DIR="${{ROOT}}/cache/chains"
TMP_DIR="${{TMPDIR:-/lustre_scratch/user_scratch/bxu2/TMPDIR/codex}}/sjcab-liftover-all.$$"
mkdir -p "${{CHAIN_DIR}}" "${{TMP_DIR}}"
trap 'rm -rf "${{TMP_DIR}}"' EXIT

case "${{REQUESTED_MODE,,}}" in
  crossmap|ucsc|auto) ;;
  *) echo "Mode must be crossmap, ucsc, or auto: ${{REQUESTED_MODE}}" >&2; exit 2 ;;
esac

ucsc_token() {{
  local value="$1" lower
  lower="${{value,,}}"
  case "${{lower}}" in
    canfam*) echo "CanFam${{value:6}}" ;;
    equcab*) echo "EquCab${{value:6}}" ;;
    bostau*) echo "BosTau${{value:6}}" ;;
    galgal*) echo "GalGal${{value:6}}" ;;
    xenla*) echo "XenLa${{value:5}}" ;;
    danrer*) echo "DanRer${{value:6}}" ;;
    susscr*) echo "SusScr${{value:6}}" ;;
    saccer*) echo "SacCer${{value:6}}" ;;
    *) echo "${{value^}}" ;;
  esac
}}

ENV_NAME=".sjcab-liftover"
activate_existing_env() {{
  if command -v micromamba >/dev/null 2>&1; then
    eval "$(micromamba shell hook -s bash)"
    micromamba activate "${{ENV_NAME}}" >/dev/null 2>&1 && return 0
  fi
  if command -v mamba >/dev/null 2>&1; then
    eval "$(mamba shell hook -s bash)"
    mamba activate "${{ENV_NAME}}" >/dev/null 2>&1 && return 0
  fi
  if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda activate "${{ENV_NAME}}" >/dev/null 2>&1 && return 0
  fi
  return 1
}}
create_and_activate_env() {{
  if command -v micromamba >/dev/null 2>&1; then
    micromamba create -y -n "${{ENV_NAME}}" -c conda-forge -c bioconda crossmap || \
      micromamba create -y -n "${{ENV_NAME}}" bioconda::ucsc-liftover
    eval "$(micromamba shell hook -s bash)"; micromamba activate "${{ENV_NAME}}"
  elif command -v mamba >/dev/null 2>&1; then
    mamba create -y -n "${{ENV_NAME}}" -c conda-forge -c bioconda crossmap || \
      mamba create -y -n "${{ENV_NAME}}" bioconda::ucsc-liftover
    eval "$(mamba shell hook -s bash)"; mamba activate "${{ENV_NAME}}"
  elif command -v conda >/dev/null 2>&1; then
    conda create -y -n "${{ENV_NAME}}" -c conda-forge -c bioconda crossmap || \
      conda create -y -n "${{ENV_NAME}}" bioconda::ucsc-liftover
    eval "$(conda shell.bash hook)"; conda activate "${{ENV_NAME}}"
  else
    echo "pixi, micromamba, mamba, or conda is required." >&2; exit 1
  fi
}}
if command -v pixi >/dev/null 2>&1; then
  PIXI_MODE=1
elif ! activate_existing_env; then
  create_and_activate_env
fi

if command -v CrossMap >/dev/null 2>&1; then
  CROSSMAP=(CrossMap)
elif command -v CrossMap.py >/dev/null 2>&1; then
  CROSSMAP=(CrossMap.py)
else
  CROSSMAP=()
fi
if command -v liftOver >/dev/null 2>&1; then LIFTOVER=(liftOver); else LIFTOVER=(); fi
if [ "${{PIXI_MODE:-0}}" = "1" ]; then LIFTOVER_MODE=pixi; fi
if [ "${{PIXI_MODE:-0}}" != 1 ] && [ "${{REQUESTED_MODE,,}}" = crossmap ] && [ "${{#CROSSMAP[@]}}" -eq 0 ]; then
  echo "CrossMap is required for mode=crossmap." >&2; exit 1
fi
if [ "${{PIXI_MODE:-0}}" != 1 ] && [ "${{REQUESTED_MODE,,}}" = ucsc ] && [ "${{#LIFTOVER[@]}}" -eq 0 ]; then
  echo "UCSC liftOver is required for mode=ucsc." >&2; exit 1
fi
if [ "${{PIXI_MODE:-0}}" != 1 ] && [ "${{REQUESTED_MODE,,}}" = auto ] && [ "${{#CROSSMAP[@]}}" -eq 0 ] && [ "${{#LIFTOVER[@]}}" -eq 0 ]; then
  echo "CrossMap or UCSC liftOver is required for mode=auto." >&2; exit 1
fi

ensure_chain() {{
  local source="$1" token name url path
  token="$(ucsc_token "${{TARGET}}")"
  name="${{source}}To${{token}}.over.chain.gz"
  path="${{CHAIN_DIR}}/${{name}}"
  url="https://hgdownload.soe.ucsc.edu/goldenPath/${{source}}/liftOver/${{name}}"
  if [ ! -s "${{path}}" ]; then curl -L "${{url}}" -o "${{path}}"; fi
  printf '%s\n' "${{path}}"
}}

run_crossmap() {{
  local chain="$1" input="$2" output="$3" fields padded mapped
  fields="$(awk 'NF {{ print NF; exit }}' "${{input}}")"
  if [ "${{fields:-0}}" -lt 12 ]; then
    padded="${{TMP_DIR}}/$(basename "${{input}}").bed12"
    mapped="${{TMP_DIR}}/$(basename "${{output}}").mapped"
    awk 'BEGIN {{ OFS="\t" }} NF {{ print $0, ".", 0, "+", $2, $3, 0, 1, $3-$2, 0 }}' "${{input}}" > "${{padded}}"
    if [ "${{LIFTOVER_MODE:-}}" = pixi ]; then
      pixi exec --channel conda-forge --channel bioconda --spec crossmap CrossMap bed "${{chain}}" "${{padded}}" "${{mapped}}"
    else
      "${{CROSSMAP[@]}}" bed "${{chain}}" "${{padded}}" "${{mapped}}"
    fi
    awk -v n="${{fields}}" 'BEGIN {{ OFS="\t" }} NF {{ for (i=1; i<=n; i++) {{ printf "%s%s", $i, (i==n ? ORS : OFS) }} }}' "${{mapped}}" > "${{output}}"
  elif [ "${{LIFTOVER_MODE:-}}" = pixi ]; then
    pixi exec --channel conda-forge --channel bioconda --spec crossmap CrossMap bed "${{chain}}" "${{input}}" "${{output}}"
  else
    "${{CROSSMAP[@]}}" bed "${{chain}}" "${{input}}" "${{output}}"
  fi
}}

run_ucsc() {{
  local chain="$1" input="$2" output="$3" fields padded mapped unmapped original
  fields="$(awk 'NF {{ print NF; exit }}' "${{input}}")"
  if [ "${{fields:-0}}" -lt 12 ]; then
    padded="${{TMP_DIR}}/$(basename "${{input}}").bed12"
    mapped="${{TMP_DIR}}/$(basename "${{output}}").mapped"
    unmapped="${{TMP_DIR}}/$(basename "${{output}}").unmap"
    original="${{TMP_DIR}}/$(basename "${{input}}").original"
    awk -v padded="${{padded}}" -v original="${{original}}" 'BEGIN {{ OFS="\t" }} NF {{
      id="sjcab_liftover_" NR
      print $1, $2, $3, id, 0, "+", $2, $3, 0, 1, $3-$2, 0 > padded
      printf "%s", id > original
      for (i=4; i<=NF; i++) printf "%s%s", OFS, $i > original
      printf "%s", ORS > original
    }}' "${{input}}"
    if [ "${{LIFTOVER_MODE:-}}" = "pixi" ]; then
      pixi exec --channel conda-forge --channel bioconda --spec ucsc-liftover liftOver "${{padded}}" "${{chain}}" "${{mapped}}" "${{unmapped}}"
    else
      "${{LIFTOVER[@]}}" "${{padded}}" "${{chain}}" "${{mapped}}" "${{unmapped}}"
    fi
    for source in "${{mapped}}" "${{unmapped}}"; do
      target="${{output}}"
      [ "${{source}}" = "${{unmapped}}" ] && target="${{output}}.unmap"
      awk -v n="${{fields}}" -v original="${{original}}" 'BEGIN {{
        OFS="\t"
        while ((getline line < original) > 0) {{
          m=split(line, a, "\t")
          for (i=2; i<=m; i++) saved[a[1], i-1]=a[i]
        }}
        close(original)
      }} /^#/ {{ print; next }} NF {{
        printf "%s\t%s\t%s", $1, $2, $3
        for (i=4; i<=n; i++) printf "%s%s", OFS, saved[$4, i]
        print ""
      }}' "${{source}}" > "${{target}}"
    done
  else
    "${{LIFTOVER[@]}}" "${{input}}" "${{chain}}" "${{output}}" "${{output}}.unmap"
  fi
}}

run_liftover() {{
  local source="$1" input="$2" output="$3" chain
  chain="$(ensure_chain "${{source}}")"
  if [ "${{LIFTOVER_MODE:-}}" = pixi ]; then
    if [ "${{REQUESTED_MODE,,}}" = ucsc ]; then
      run_ucsc "${{chain}}" "${{input}}" "${{output}}"
    else
      run_crossmap "${{chain}}" "${{input}}" "${{output}}"
    fi
  elif [ "${{REQUESTED_MODE,,}}" = crossmap ] || {{ [ "${{REQUESTED_MODE,,}}" = auto ] && [ "${{#CROSSMAP[@]}}" -gt 0 ]; }}; then
    run_crossmap "${{chain}}" "${{input}}" "${{output}}"
  else
    run_ucsc "${{chain}}" "${{input}}" "${{output}}"
  fi
}}

process_file() {{
  local source="$1" input="$2" output="$3" compressed=0 plain
  if [[ "${{output}}" == *.gz ]]; then compressed=1; plain="${{output%.gz}}"; else plain="${{output}}"; fi
  if [ -s "${{output}}" ] || [ -s "${{plain}}" ]; then return 0; fi
  mkdir -p "$(dirname "${{plain}}")"
  if [[ "${{input}}" == *.gz ]]; then
    local unpacked="${{TMP_DIR}}/$(basename "${{plain}}")"
    gzip -dc "${{input}}" > "${{unpacked}}"
    run_liftover "${{source}}" "${{unpacked}}" "${{plain}}"
  else
    run_liftover "${{source}}" "${{input}}" "${{plain}}"
  fi
  [ -e "${{plain}}.unmap" ] || : > "${{plain}}.unmap"
  if [ "${{compressed}}" -eq 1 ]; then gzip -f "${{plain}}"; fi
}}

for input in "${{ROOT}}"/blacklists/hg38*.bed "${{ROOT}}"/cgi/hg38*.bed; do
  [ -f "${{input}}" ] || continue
  directory="$(dirname "${{input}}")"
  base="$(basename "${{input}}")"
  output_base="${{base#hg38}}"
  process_file hg38 "${{input}}" "${{directory}}/${{TARGET}}${{output_base}}"
done

for input in "${{ROOT}}"/chromhmm/hg38/**/*.bed "${{ROOT}}"/chromhmm/hg38/**/*.bed.gz; do
  [ -f "${{input}}" ] || continue
  relative="${{input#${{ROOT}}/chromhmm/hg38/}}"
  process_file hg38 "${{input}}" "${{ROOT}}/chromhmm/${{TARGET}}/${{relative}}"
done

for input in "${{ROOT}}"/segway/hg19/*.bed.gz "${{ROOT}}"/segway/hg19/interpreted/*.bed.gz; do
  [ -f "${{input}}" ] || continue
  relative="$(basename "${{input}}" .bed.gz)"
  process_file hg19 "${{input}}" "${{ROOT}}/segway/${{TARGET}}/${{relative}}.${{TARGET}}lift.bed.gz"
done
'''.format(
        root=_shell_quote(str(root)),
        target=_shell_quote(_genome(target_genome)),
        mode=selected_mode,
    )
    temporary = output.with_name(output.name + ".tmp")
    try:
        temporary.write_text(script, encoding="utf-8")
        temporary.chmod(0o755)
        temporary.replace(output)
    finally:
        if temporary.exists():
            temporary.unlink()
    return output


def _genome(value: str) -> str:
    genome = str(value).strip()
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]*", genome):
        raise ValueError("Genome must be a simple UCSC-style build name.")
    return genome.lower()


def _ucsc_token(genome: str) -> str:
    """Return the canonical mixed-case token used by UCSC chain files."""

    prefixes = {
        "canf": "CanF",
        "equcab": "EquCab",
        "bostau": "BosTau",
        "galgal": "GalGal",
        "xenla": "XenLa",
        "danrer": "DanRer",
        "susscr": "SusScr",
        "saccer": "SacCer",
        "orycun": "OryCun",
        "mondom": "MonDom",
        "pantro": "PanTro",
        "gorgor": "GorGor",
        "ponabe": "PonAbe",
        "caljac": "CalJac",
        "chlsab": "ChlSab",
        "tetnig": "TetNig",
    }
    lowered = genome.lower()
    for prefix, canonical in prefixes.items():
        if lowered.startswith(prefix):
            return canonical + genome[len(prefix) :]
    return genome[:1].upper() + genome[1:]


def _shell_quote(value: str) -> str:
    return "'{}'".format(value.replace("'", "'\\''"))
