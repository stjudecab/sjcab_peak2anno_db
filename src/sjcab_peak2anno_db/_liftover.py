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
    env_name: str = "sjcab-liftover",
    install_root: Optional[PathLike] = None,
    rename_genome_prefix: bool = False,
    cache_dir: Optional[PathLike] = None,
) -> str:
    """Return a bash script that lifts BED and BED.GZ files in a directory."""

    source = _genome(source_genome)
    target = _genome(target_genome)
    if source == target:
        raise ValueError("Source and target genomes must differ for liftover.")
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
    return """#!/usr/bin/env bash
set -euo pipefail

ENV_NAME={env_name}
TARGET={target}
INPUT_DIR={input_dir}
OUTPUT_DIR={output_dir}
INSTALL_ROOT={install_root}
CHAIN_DIR={chain_dir}
CHAIN="${{CHAIN_DIR}}/{chain_name}"
TMP_DIR="${{TMPDIR:-/lustre_scratch/user_scratch/bxu2/TMPDIR/codex}}/sjcab-liftover.$$"
mkdir -p "${{TMP_DIR}}"
trap 'rm -rf "${{TMP_DIR}}"' EXIT

mkdir -p "${{CHAIN_DIR}}" "${{OUTPUT_DIR}}"
if [ ! -s "${{CHAIN}}" ]; then
  curl -L "{chain_url}" -o "${{CHAIN}}"
fi

activate_existing_env() {{
  if command -v micromamba >/dev/null 2>&1; then
    hook="$(micromamba shell hook -s bash 2>/dev/null)" && eval "${{hook}}"
    micromamba activate "${{ENV_NAME}}" >/dev/null 2>&1 && return 0
  fi
  if command -v conda >/dev/null 2>&1; then
    hook="$(conda shell.bash hook 2>/dev/null)" && eval "${{hook}}"
    conda activate "${{ENV_NAME}}" >/dev/null 2>&1 && return 0
  fi
  if command -v mamba >/dev/null 2>&1; then
    hook="$(mamba shell hook -s bash 2>/dev/null)" && eval "${{hook}}"
    mamba activate "${{ENV_NAME}}" >/dev/null 2>&1 && return 0
  fi
  return 1
}}

create_and_activate_env() {{
  if command -v micromamba >/dev/null 2>&1; then
    if micromamba create -y -n "${{ENV_NAME}}" bioconda::ucsc-liftover; then
      eval "$(micromamba shell hook -s bash)"
      micromamba activate "${{ENV_NAME}}"
      return 0
    fi
    micromamba create -y -n "${{ENV_NAME}}" -c conda-forge -c bioconda crossmap
    eval "$(micromamba shell hook -s bash)"
    micromamba activate "${{ENV_NAME}}"
  elif command -v mamba >/dev/null 2>&1; then
    if mamba create -y -n "${{ENV_NAME}}" bioconda::ucsc-liftover; then
      eval "$(mamba shell hook -s bash)"
      mamba activate "${{ENV_NAME}}"
      return 0
    fi
    mamba create -y -n "${{ENV_NAME}}" -c conda-forge -c bioconda crossmap
    eval "$(mamba shell hook -s bash)"
    mamba activate "${{ENV_NAME}}"
  elif command -v conda >/dev/null 2>&1; then
    if conda create -y -n "${{ENV_NAME}}" bioconda::ucsc-liftover; then
      eval "$(conda shell.bash hook)"
      conda activate "${{ENV_NAME}}"
      return 0
    fi
    conda create -y -n "${{ENV_NAME}}" -c conda-forge -c bioconda crossmap
    eval "$(conda shell.bash hook)"
    conda activate "${{ENV_NAME}}"
  else
    echo "micromamba, mamba, or conda is required to install liftOver or CrossMap" >&2
    exit 1
  fi
}}

if ! activate_existing_env; then
  create_and_activate_env
fi

if command -v liftOver >/dev/null 2>&1; then
  LIFTOVER=(liftOver)
  LIFTOVER_MODE=ucsc
elif command -v CrossMap >/dev/null 2>&1; then
  CROSSMAP=(CrossMap)
  LIFTOVER_MODE=crossmap
elif command -v CrossMap.py >/dev/null 2>&1; then
  CROSSMAP=(CrossMap.py)
  LIFTOVER_MODE=crossmap
else
  echo "liftOver and CrossMap are not available in conda environment ${{ENV_NAME}}." >&2
  exit 1
fi

run_liftover() {{
  if [ "${{LIFTOVER_MODE}}" = "ucsc" ]; then
    "${{LIFTOVER[@]}}" "$@"
  else
    "${{CROSSMAP[@]}}" bed "$@"
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
        chain_url=chain_url,
        chain_dir=_shell_quote(str(chain_cache_dir)),
        env_name=_shell_quote(env_name),
        input_dir=_shell_quote(str(Path(input_dir).expanduser())),
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
    )


def write_liftover_script(
    input_dir: PathLike,
    output_dir: PathLike,
    source_genome: str,
    target_genome: str,
    script_path: PathLike,
    env_name: str = "sjcab-liftover",
    install_root: Optional[PathLike] = None,
    rename_genome_prefix: bool = False,
    cache_dir: Optional[PathLike] = None,
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
            ),
            encoding="utf-8",
        )
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
