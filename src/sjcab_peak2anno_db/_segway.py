"""Download Segway encyclopedia hg19 BED resources."""

from __future__ import annotations

import json
import re
import urllib.parse
import urllib.request
from dataclasses import dataclass
from difflib import SequenceMatcher
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

from ._download import ProgressCallback, download_file, report_progress
from ._download_log import record_download_url
from ._registry import UnknownResourceError, user_data_dir

PathLike = Union[str, Path]

SEGWAY_DIR_NAME = "segway"
SEGWAY_GENOMES = ("hg19",)
SEGWAY_URL = (
    "https://www.encodeproject.org/publications/"
    "94941f71-80c8-43d2-809b-25161efc3be0/"
)
SEGWAY_FALLBACK_URL = "https://noble.gs.washington.edu/proj/encyclopedia/"
SEGWAY_INTERPRETED_DIR = "interpreted/"
SEGWAY_CHAIN_ROOT_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/"
SEGWAY_METADATA_TSV_NAME = "metadata.tsv"
SEGWAY_METADATA_FIELDS = "encodeID,biosample_term_name,biosample_type,biosample_term_id,filename,url,encode_url".split(",")
SEGWAY_DOWNLOAD_MANIFEST_NAME = "segway_downloads.tsv"
SEGWAY_AGGREGATE_FILES = (
    "segway_encyclopedia.bed.gz",
    "caas.bed.gz",
    "label_info.tab",
)
SEGWAY_DEFAULT_FILE = "segway_encyclopedia.bed.gz"
ENCODE_BASE_URL = "https://www.encodeproject.org"
ENCODE_FILE_ACCESSION_RE = re.compile(r"^ENCFF[0-9]{3}[A-Z]{3}$")
ENCODE_GENERIC_SEGWAY_BASENAMES = {
    "ann.bed",
    "ann.bed.gz",
    "segway.bed",
    "segway.bed.gz",
    "recolored.bed",
    "recolored.bed.gz",
}
_BUNDLED_SEGWAY_FILES = None  # type: Optional[Dict[str, Tuple[SegwayFile, ...]]]


@dataclass(frozen=True)
class SegwayFile:
    """Segway downloadable file metadata."""

    name: str
    url: str
    encode_url: str = ""
    kind: str

    @property
    def key(self) -> str:
        if self.name in SEGWAY_AGGREGATE_FILES:
            return _strip_known_suffix(self.name)
        return _strip_known_suffix(self.name)

    @property
    def relative_path(self) -> Path:
        return Path(self.name)

    @property
    def search_text(self) -> str:
        return " ".join((self.name, self.key, self.kind))


def segway_root(data_dir: Optional[PathLike] = None) -> Path:
    """Return the local Segway cache root."""

    return user_data_dir(data_dir) / SEGWAY_DIR_NAME


def segway_dir(
    data_dir: Optional[PathLike] = None,
    genome: str = "hg19",
) -> Path:
    """Return the local cache directory for Segway resources."""

    genome = _normalize_genome(genome)
    return segway_root(data_dir) / genome


def segway_filename(name: str = SEGWAY_DEFAULT_FILE) -> str:
    """Return a normalized Segway resource filename."""

    value = name.strip()
    if not value:
        raise UnknownResourceError("Segway file names cannot be empty.")
    if value in SEGWAY_AGGREGATE_FILES:
        return value
    if "/" in value or "\\" in value:
        raise UnknownResourceError("Segway file names should not contain directories.")
    if value.endswith(".bed.gz") or value.endswith(".tab"):
        return value
    return "{}.bed.gz".format(value)


def segway_path(
    name: str = SEGWAY_DEFAULT_FILE,
    data_dir: Optional[PathLike] = None,
    genome: str = "hg19",
) -> Path:
    """Return the local cache path for one Segway resource."""

    genome = _normalize_genome(genome)
    filename = segway_filename(name)
    return segway_dir(data_dir, genome) / filename


def segway_metadata_tsv_path(
    data_dir: Optional[PathLike] = None,
    genome: str = "hg19",
) -> Path:
    """Return the local Segway metadata TSV path."""

    _normalize_genome(genome)
    return segway_root(data_dir) / SEGWAY_METADATA_TSV_NAME


def segway_download_manifest_path(data_dir: Optional[PathLike] = None) -> Path:
    """Return the installed Segway download manifest TSV path."""

    return segway_root(data_dir) / SEGWAY_DOWNLOAD_MANIFEST_NAME


def segway_files(
    source_url: str = SEGWAY_URL,
    fallback_url: str = SEGWAY_FALLBACK_URL,
) -> Tuple[SegwayFile, ...]:
    """Return Segway files discovered from the primary or fallback source page."""

    entries, _source_page = _discover_files(source_url, fallback_url)
    return entries


def segway_available_names(
    source_url: str = SEGWAY_URL,
    fallback_url: str = SEGWAY_FALLBACK_URL,
) -> Tuple[str, ...]:
    """Return available Segway cell-type/sample names."""

    return tuple(
        entry.key
        for entry in segway_files(source_url, fallback_url)
        if entry.kind == "celltype"
    )


def segway_url(
    name: str = SEGWAY_DEFAULT_FILE,
    source_url: str = SEGWAY_URL,
    fallback_url: str = SEGWAY_FALLBACK_URL,
) -> str:
    """Return the URL for one Segway resource."""

    filename = segway_filename(name)
    entries = segway_files(source_url, fallback_url)
    aliases = _entry_aliases(entries)
    entry = aliases.get(_normalize_text(_strip_known_suffix(name)))
    if entry is not None:
        return entry.encode_url or entry.url
    for entry in entries:
        if entry.name == filename:
            return entry.encode_url or entry.url
    raise UnknownResourceError("No Segway file matched {!r}.".format(name))


def download_segway(
    data_dir: Optional[PathLike] = None,
    genome: str = "hg19",
    names: Optional[Union[str, Iterable[str]]] = None,
    tissue: Optional[Union[str, Iterable[str]]] = None,
    cellline: Optional[Union[str, Iterable[str]]] = None,
    all_celltypes: bool = False,
    include_encyclopedia: bool = False,
    include_caas: bool = False,
    include_label_info: bool = False,
    overwrite: bool = False,
    source_url: str = SEGWAY_URL,
    fallback_url: str = SEGWAY_FALLBACK_URL,
    progress: Optional[ProgressCallback] = None,
) -> Mapping[str, Path]:
    """Download Segway encyclopedia resources.

    Segway encyclopedia downloads are hg19 only. With no selector, the
    cell type-agnostic ``segway_encyclopedia.bed.gz`` file is downloaded.
    Select cell-type BED files by sample names such as ``"GM12878"`` or fuzzy
    tissue/cell-line keywords. Set ``all_celltypes`` to download every
    discovered cell-type annotation.
    """

    genome = _normalize_genome(genome)
    selected, source_page = _select_files_from_sources(
        source_url,
        fallback_url,
        progress=progress,
        names=names,
        tissue=tissue,
        cellline=cellline,
        all_celltypes=all_celltypes,
        include_encyclopedia=include_encyclopedia,
        include_caas=include_caas,
        include_label_info=include_label_info,
    )

    target_dir = segway_dir(data_dir, genome)
    target_dir.mkdir(parents=True, exist_ok=True)

    outputs = {}  # type: Dict[str, Path]
    metadata_entries = {}  # type: Dict[str, SegwayFile]
    for entry in selected:
        destination = target_dir / entry.relative_path
        destination.parent.mkdir(parents=True, exist_ok=True)
        if overwrite or not destination.exists():
            if progress is None:
                _download_file(entry.encode_url or entry.url, destination)
            else:
                _download_file(entry.encode_url or entry.url, destination, progress=progress)
            record_download_url(entry.encode_url or entry.url, data_dir=data_dir, destination=destination)
        outputs[entry.key] = destination
        metadata_entries[entry.key] = entry

    _write_segway_metadata_tsv(
        outputs,
        metadata_entries,
        genome=genome,
        data_dir=data_dir,
        progress=progress,
    )
    _install_segway_download_manifest(data_dir)
    record_download_url(source_page, data_dir=data_dir)
    return outputs


def install_segway(
    data_dir: Optional[PathLike] = None,
    genome: str = "hg19",
    names: Optional[Union[str, Iterable[str]]] = None,
    tissue: Optional[Union[str, Iterable[str]]] = None,
    cellline: Optional[Union[str, Iterable[str]]] = None,
    all_celltypes: bool = False,
    include_encyclopedia: bool = False,
    include_caas: bool = False,
    include_label_info: bool = False,
    overwrite: bool = False,
    source_url: str = SEGWAY_URL,
    fallback_url: str = SEGWAY_FALLBACK_URL,
    progress: Optional[ProgressCallback] = None,
) -> Mapping[str, Path]:
    """Install Segway resources into the configured user data directory."""

    return download_segway(
        data_dir=data_dir,
        genome=genome,
        names=names,
        tissue=tissue,
        cellline=cellline,
        all_celltypes=all_celltypes,
        include_encyclopedia=include_encyclopedia,
        include_caas=include_caas,
        include_label_info=include_label_info,
        overwrite=overwrite,
        source_url=source_url,
        fallback_url=fallback_url,
        progress=progress,
    )


def segway_liftover_script(
    input_dir: PathLike,
    output_dir: PathLike,
    target_genome: str,
    env_name: str = "segway-liftover",
    install_root: Optional[PathLike] = None,
) -> str:
    """Return a bash script that lifts local hg19 Segway BEDs to another genome."""

    target = _normalize_target_genome(target_genome)
    chain_name = "hg19To{}.over.chain.gz".format(_ucsc_target_token(target))
    chain_url = urllib.parse.urljoin(SEGWAY_CHAIN_ROOT_URL, chain_name)
    install_root_text = (
        _shell_quote(str(Path(install_root).expanduser()))
        if install_root is not None
        else "''"
    )
    return """#!/usr/bin/env bash
set -euo pipefail

ENV_NAME={env_name}
INPUT_DIR={input_dir}
OUTPUT_DIR={output_dir}
INSTALL_ROOT={install_root}
CHAIN_DIR="${{OUTPUT_DIR}}/chains"
CHAIN="${{CHAIN_DIR}}/{chain_name}"

mkdir -p "${{CHAIN_DIR}}" "${{OUTPUT_DIR}}"
curl -L "{chain_url}" -o "${{CHAIN}}"

activate_existing_env() {{
  if command -v micromamba >/dev/null 2>&1; then
    hook="$(micromamba shell hook -s bash 2>/dev/null)" && eval "${{hook}}"
    if micromamba activate "${{ENV_NAME}}" >/dev/null 2>&1; then
      echo "Activated existing micromamba environment: ${{ENV_NAME}}" >&2
      return 0
    fi
  fi

  if command -v conda >/dev/null 2>&1; then
    hook="$(conda shell.bash hook 2>/dev/null)" && eval "${{hook}}"
    if conda activate "${{ENV_NAME}}" >/dev/null 2>&1; then
      echo "Activated existing conda environment: ${{ENV_NAME}}" >&2
      return 0
    fi
  fi

  if command -v mamba >/dev/null 2>&1; then
    hook="$(mamba shell hook -s bash 2>/dev/null)" && eval "${{hook}}"
    if mamba activate "${{ENV_NAME}}" >/dev/null 2>&1; then
      echo "Activated existing mamba environment: ${{ENV_NAME}}" >&2
      return 0
    fi
  fi

  return 1
}}

create_and_activate_env() {{
  if command -v micromamba >/dev/null 2>&1; then
    micromamba create -y -n "${{ENV_NAME}}" -c conda-forge -c bioconda crossmap
    eval "$(micromamba shell hook -s bash)"
    micromamba activate "${{ENV_NAME}}"
  elif command -v mamba >/dev/null 2>&1; then
    mamba create -y -n "${{ENV_NAME}}" -c conda-forge -c bioconda crossmap
    if command -v conda >/dev/null 2>&1; then
      eval "$(conda shell.bash hook)"
      conda activate "${{ENV_NAME}}"
    else
      eval "$(mamba shell hook -s bash)"
      mamba activate "${{ENV_NAME}}"
    fi
  elif command -v conda >/dev/null 2>&1; then
    conda create -y -n "${{ENV_NAME}}" -c conda-forge -c bioconda crossmap
    eval "$(conda shell.bash hook)"
    conda activate "${{ENV_NAME}}"
  else
    echo "micromamba, mamba, or conda is required to install CrossMap" >&2
    exit 1
  fi
}}

if ! activate_existing_env; then
  create_and_activate_env
fi

if command -v CrossMap >/dev/null 2>&1; then
  CROSSMAP=(CrossMap)
elif command -v CrossMap.py >/dev/null 2>&1; then
  CROSSMAP=(CrossMap.py)
else
  echo "CrossMap is not available in conda environment ${{ENV_NAME}}." >&2
  echo "Install it there manually or remove the env so this script can recreate it." >&2
  exit 1
fi

for bed in "${{INPUT_DIR}}"/*.bed.gz "${{INPUT_DIR}}"/interpreted/*.bed.gz; do
  [ -e "${{bed}}" ] || continue
  rel="$(basename "${{bed}}")"
  out="${{OUTPUT_DIR}}/${{rel%.bed.gz}}.{target}lift.bed"
  if [ -s "${{out}}.gz" ]; then
    echo "Skipping existing lifted file: ${{out}}.gz" >&2
    continue
  fi
  mkdir -p "$(dirname "${{out}}")"
  "${{CROSSMAP[@]}}" bed "${{CHAIN}}" "${{bed}}" "${{out}}"
  gzip -f "${{out}}"
done

if [ -n "${{INSTALL_ROOT}}" ]; then
  INSTALL_DIR="${{INSTALL_ROOT}}/{target}"
  mkdir -p "${{INSTALL_DIR}}"
  for lifted in "${{OUTPUT_DIR}}"/*.bed.gz; do
    [ -e "${{lifted}}" ] || continue
    rel="$(basename "${{lifted}}")"
    cp -p "${{lifted}}" "${{INSTALL_DIR}}/${{rel}}"
  done
  echo "Lifted Segway files are organized under ${{INSTALL_DIR}}" >&2
else
  echo "Lifted Segway files are organized under ${{OUTPUT_DIR}}" >&2
fi
""".format(
        chain_name=chain_name,
        chain_url=chain_url,
        env_name=_shell_quote(env_name),
        input_dir=_shell_quote(str(Path(input_dir).expanduser())),
        install_root=install_root_text,
        output_dir=_shell_quote(str(Path(output_dir).expanduser())),
        target=target,
    )


def write_segway_liftover_script(
    data_dir: Optional[PathLike] = None,
    target_genome: str = "hg38",
    script_path: Optional[PathLike] = None,
    install_root: Optional[PathLike] = None,
) -> Path:
    """Write a Segway hg19 liftover helper script and return its path."""

    target = _normalize_target_genome(target_genome)
    root = segway_root(data_dir)
    input_dir = root / "hg19"
    if install_root is None:
        output_dir = root / target
    else:
        output_dir = root / ".liftover" / target
    if script_path is None:
        script_path = root / "liftover_hg19_to_{}.sh".format(target)
    script = segway_liftover_script(
        input_dir,
        output_dir,
        target,
        install_root=install_root,
    )
    output = Path(script_path).expanduser()
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output.with_name(output.name + ".tmp")
    try:
        with tmp_path.open("w", encoding="utf-8") as handle:
            handle.write(script)
        tmp_path.chmod(0o755)
        tmp_path.replace(output)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()
    return output


def _discover_files(
    source_url: str,
    fallback_url: str,
    progress: Optional[ProgressCallback] = None,
) -> Tuple[Tuple[SegwayFile, ...], str]:
    errors = []  # type: List[str]
    source_urls = _unique_source_urls(source_url, fallback_url)
    if _all_sources_are_builtin(source_urls):
        entries = _combined_bundled_manifest_entries(source_urls)
        if entries:
            return entries, source_urls[0]

    for index, page_url in enumerate(source_urls):
        try:
            entries, source_page = _discover_files_for_source(page_url)
        except Exception as exc:
            reason = "{} ({})".format(page_url, exc)
            errors.append(reason)
            if index == 0 and fallback_url:
                report_progress(
                    progress,
                    "download-segway: primary source unavailable: {}; "
                    "trying fallback {}".format(reason, fallback_url),
                )
            continue
        if entries:
            if index > 0:
                report_progress(
                    progress,
                    "download-segway: using fallback source {}".format(page_url),
                )
            return entries, source_page
        reason = "{} (no Segway download links found)".format(page_url)
        errors.append(reason)
        if index == 0 and fallback_url:
            report_progress(
                progress,
                "download-segway: primary source unavailable: {}; "
                "trying fallback {}".format(reason, fallback_url),
            )
    raise UnknownResourceError(
        "Could not discover Segway downloads from {} or {}. {}".format(
            source_url,
            fallback_url,
            "; ".join(errors),
        )
    )


def _select_files_from_sources(
    source_url: str,
    fallback_url: str,
    progress: Optional[ProgressCallback],
    names: Optional[Union[str, Iterable[str]]],
    tissue: Optional[Union[str, Iterable[str]]],
    cellline: Optional[Union[str, Iterable[str]]],
    all_celltypes: bool,
    include_encyclopedia: bool,
    include_caas: bool,
    include_label_info: bool,
) -> Tuple[Tuple[SegwayFile, ...], str]:
    errors = []  # type: List[str]
    source_urls = _unique_source_urls(source_url, fallback_url)
    if _all_sources_are_builtin(source_urls):
        entries = _combined_bundled_manifest_entries(source_urls)
        if entries:
            try:
                selected = _select_files(
                    entries,
                    names=names,
                    tissue=tissue,
                    cellline=cellline,
                    all_celltypes=all_celltypes,
                    include_encyclopedia=include_encyclopedia,
                    include_caas=include_caas,
                    include_label_info=include_label_info,
                )
            except UnknownResourceError as exc:
                raise UnknownResourceError(
                    "No Segway files matched the request from {}. {}".format(
                        " or ".join(source_urls),
                        exc,
                    )
                )
            if selected:
                fallback_source = _selected_fallback_source_url(selected, source_urls)
                selected_source = source_urls[0]
                if fallback_source is not None:
                    report_progress(
                        progress,
                        "download-segway: using fallback source {}".format(
                            fallback_source
                        ),
                    )
                    fallback_source_name = _manifest_source_name(fallback_source)
                    if all(
                        _manifest_source_name(entry.url) == fallback_source_name
                        for entry in selected
                    ):
                        selected_source = fallback_source
                return selected, selected_source

    for index, page_url in enumerate(source_urls):
        try:
            entries, source_page = _discover_files_for_source(page_url)
        except Exception as exc:
            reason = "{} ({})".format(page_url, exc)
            errors.append(reason)
            if index == 0 and len(source_urls) > 1:
                report_progress(
                    progress,
                    "download-segway: primary source unavailable: {}; "
                    "trying fallback {}".format(reason, fallback_url),
                )
            continue

        try:
            selected = _select_files(
                entries,
                names=names,
                tissue=tissue,
                cellline=cellline,
                all_celltypes=all_celltypes,
                include_encyclopedia=include_encyclopedia,
                include_caas=include_caas,
                include_label_info=include_label_info,
            )
        except UnknownResourceError as exc:
            reason = "{} ({})".format(page_url, exc)
            errors.append(reason)
            if (
                index == 0
                and len(source_urls) > 1
                and not _is_builtin_encode_source(page_url)
            ):
                report_progress(
                    progress,
                    "download-segway: primary source had no matching Segway files; "
                    "trying fallback {}".format(fallback_url),
                )
            continue

        if selected:
            if index > 0:
                report_progress(
                    progress,
                    "download-segway: using fallback source {}".format(page_url),
                )
            return selected, source_page

        reason = "{} (no Segway files matched the request)".format(page_url)
        errors.append(reason)
        if (
            index == 0
            and len(source_urls) > 1
            and not _is_builtin_encode_source(page_url)
        ):
            report_progress(
                progress,
                "download-segway: primary source had no matching Segway files; "
                "trying fallback {}".format(fallback_url),
            )

    raise UnknownResourceError(
        "No Segway files matched the request from {}. {}".format(
            " or ".join(source_urls),
            "; ".join(errors),
        )
    )


def _discover_files_for_source(page_url: str) -> Tuple[Tuple[SegwayFile, ...], str]:
    entries = _bundled_manifest_entries_for_url(page_url)
    if entries:
        return entries, page_url
    if _is_builtin_source_url(page_url):
        return tuple(), page_url
    return _discover_files_from_page(page_url), page_url


def _combined_bundled_manifest_entries(
    source_urls: Sequence[str],
) -> Tuple[SegwayFile, ...]:
    entries = []  # type: List[SegwayFile]
    seen = set()
    for page_url in source_urls:
        for entry in _bundled_manifest_entries_for_url(page_url):
            key = (entry.name, entry.kind)
            if key in seen:
                continue
            entries.append(entry)
            seen.add(key)
    return tuple(entries)


def _bundled_manifest_entries_for_url(page_url: str) -> Tuple[SegwayFile, ...]:
    if not _is_builtin_source_url(page_url):
        return tuple()
    source_name = _manifest_source_name(page_url)
    if source_name is None:
        return tuple()
    return _bundled_manifest_entries().get(source_name, tuple())


def _bundled_manifest_entries() -> Dict[str, Tuple[SegwayFile, ...]]:
    global _BUNDLED_SEGWAY_FILES
    if _BUNDLED_SEGWAY_FILES is not None:
        return _BUNDLED_SEGWAY_FILES

    path = _bundled_segway_manifest_path()
    entries = {}  # type: Dict[str, List[SegwayFile]]
    if path.exists():
        with path.open("r", encoding="utf-8") as handle:
            header = handle.readline().rstrip("\n").split("\t")
            if header == ["source", "name", "kind", "url", "encode_url"]:
                for line in handle:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) != 5:
                        continue
                    source, name, kind, url, encode_url = fields
                    entries.setdefault(source, []).append(
                        SegwayFile(name=name, kind=kind, url=url, encode_url=encode_url)
                    )
            elif header == ["encodeID", "biosample_term_name", "biosample_type", "biosample_term_id", "filename", "url", "encode_url"]:
                for line in handle:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) != 7:
                        continue
                    encode_id, biosample_term_name, biosample_type, biosample_term_id, filename, url, encode_url = fields
                    entries.setdefault("encode", []).append(
                        SegwayFile(name=filename, kind="celltype", url=url, encode_url=encode_url)
                    )

    _BUNDLED_SEGWAY_FILES = {
        source: tuple(values) for source, values in entries.items()
    }
    return _BUNDLED_SEGWAY_FILES


def _manifest_source_name(page_url: str) -> Optional[str]:
    parsed = urllib.parse.urlparse(page_url)
    netloc = parsed.netloc.lower()
    if netloc.endswith("encodeproject.org"):
        return "encode"
    if netloc == "noble.gs.washington.edu":
        return "washington"
    return None


def _is_builtin_source_url(page_url: str) -> bool:
    return _normalized_source_url(page_url) in {
        _normalized_source_url(SEGWAY_URL),
        _normalized_source_url(SEGWAY_FALLBACK_URL),
        _normalized_source_url(ENCODE_BASE_URL),
    }


def _is_builtin_encode_source(page_url: str) -> bool:
    return (
        _manifest_source_name(page_url) == "encode"
        and _is_builtin_source_url(page_url)
    )


def _all_sources_are_builtin(source_urls: Sequence[str]) -> bool:
    return bool(source_urls) and all(_is_builtin_source_url(url) for url in source_urls)


def _selected_fallback_source_url(
    selected: Sequence[SegwayFile],
    source_urls: Sequence[str],
) -> Optional[str]:
    if not source_urls:
        return None
    primary_source = _manifest_source_name(source_urls[0])
    for page_url in source_urls[1:]:
        source_name = _manifest_source_name(page_url)
        if source_name is None or source_name == primary_source:
            continue
        if any(_manifest_source_name(entry.url) == source_name for entry in selected):
            return page_url
    return None


def _unique_source_urls(source_url: str, fallback_url: str) -> Tuple[str, ...]:
    urls = []  # type: List[str]
    for url in (source_url, fallback_url):
        if url and url not in urls:
            urls.append(url)
    return tuple(urls)


def _normalized_source_url(url: str) -> str:
    parsed = urllib.parse.urlparse(url)
    path = _ensure_trailing_slash(parsed.path)
    return urllib.parse.urlunparse(
        parsed._replace(path=path, params="", query="", fragment="")
    )


def _discover_files_from_page(page_url: str) -> Tuple[SegwayFile, ...]:
    if _is_encode_url(page_url):
        return _discover_files_from_encode(page_url)

    html = _fetch_text(page_url)
    if _is_human_verification_page(html):
        raise ValueError("human verification/captcha page returned")
    entries = []  # type: List[SegwayFile]
    seen = set()

    for href in _html_links(html):
        filename = urllib.parse.unquote(Path(urllib.parse.urlparse(href).path).name)
        if filename in SEGWAY_AGGREGATE_FILES:
            _append_entry(
                entries,
                seen,
                filename,
                urllib.parse.urljoin(page_url, href),
                encode_url=_aggregate_kind(filename),
            )

    for directory_url in _interpreted_directory_urls(html, page_url):
        try:
            directory_html = _fetch_text(directory_url)
        except Exception:
            continue
        for href in _html_links(directory_html):
            filename = urllib.parse.unquote(Path(urllib.parse.urlparse(href).path).name)
            if filename.endswith(".bed.gz") and filename not in SEGWAY_AGGREGATE_FILES:
                _append_entry(
                    entries,
                    seen,
                    filename,
                    urllib.parse.urljoin(directory_url, href),
                    "celltype",
                )

    for href in _html_links(html):
        filename = urllib.parse.unquote(Path(urllib.parse.urlparse(href).path).name)
        if filename.endswith(".bed.gz") and filename not in SEGWAY_AGGREGATE_FILES:
            _append_entry(
                entries,
                seen,
                filename,
                urllib.parse.urljoin(page_url, href),
                "celltype",
            )

    return tuple(entries)


def _discover_files_from_encode(page_url: str) -> Tuple[SegwayFile, ...]:
    payload = _fetch_encode_json(_encode_json_url(page_url))
    entries = list(_encode_entries_from_document(payload, page_url))

    if not entries:
        for search_url in _encode_related_search_urls(payload, page_url):
            try:
                search_payload = _fetch_encode_json(search_url)
            except Exception:
                continue
            entries.extend(_encode_entries_from_document(search_payload, page_url))

    return _dedupe_entries(entries)


def _encode_entries_from_document(
    payload: Any,
    page_url: str,
) -> Tuple[SegwayFile, ...]:
    entries = []  # type: List[SegwayFile]
    for record in _iter_encode_file_objects(payload):
        entry = _encode_entry_from_file(record, page_url)
        if entry is not None:
            entries.append(entry)
    return _dedupe_entries(entries)


def _encode_entry_from_file(
    record: Mapping[str, Any],
    page_url: str,
) -> Optional[SegwayFile]:
    accession = _encode_file_accession(record)
    if accession is None:
        return None
    if record.get("no_file_available"):
        return None
    status = str(record.get("status", "")).lower()
    if status in {"deleted", "revoked"}:
        return None

    basename = _encode_submitted_basename(record)
    kind = _encode_kind(record, basename)
    name = _encode_local_filename(record, accession, basename, kind)
    if name is None:
        return None
    if not _is_encode_hg19_record(record, kind):
        return None
    if not _is_encode_segway_file(record, name, kind):
        return None

    return SegwayFile(
        name=name,
        url=_encode_download_url(record, accession, page_url),
        kind=kind,
    )


def _iter_encode_file_objects(value: Any) -> Iterable[Mapping[str, Any]]:
    if isinstance(value, Mapping):
        if _looks_like_encode_file(value):
            yield value
        for child in value.values():
            yield from _iter_encode_file_objects(child)
    elif isinstance(value, list):
        for child in value:
            yield from _iter_encode_file_objects(child)


def _looks_like_encode_file(value: Mapping[str, Any]) -> bool:
    types = value.get("@type", [])
    if isinstance(types, str):
        types = [types]
    if "File" in types:
        return True
    accession = value.get("accession")
    if isinstance(accession, str) and ENCODE_FILE_ACCESSION_RE.match(accession):
        return True
    href = value.get("href")
    return isinstance(href, str) and "/files/ENCFF" in href and "/@@download/" in href


def _encode_file_accession(record: Mapping[str, Any]) -> Optional[str]:
    accession = record.get("accession")
    if isinstance(accession, str) and ENCODE_FILE_ACCESSION_RE.match(accession):
        return accession
    href = record.get("href")
    if isinstance(href, str):
        match = re.search(r"/files/(ENCFF[0-9]{3}[A-Z]{3})/", href)
        if match:
            return match.group(1)
    file_id = record.get("@id")
    if isinstance(file_id, str):
        match = re.search(r"/files/(ENCFF[0-9]{3}[A-Z]{3})/", file_id)
        if match:
            return match.group(1)
    return None


def _encode_kind(record: Mapping[str, Any], basename: str) -> str:
    text = _encode_record_text(record)
    if basename in SEGWAY_AGGREGATE_FILES:
        return _aggregate_kind(basename)
    if "label_info" in text or "label info" in text:
        return "label_info"
    if "caas" in text or "conservation-associated activity" in text:
        return "caas"
    if "encyclopedia" in text and "segway" in text:
        return "encyclopedia"
    return "celltype"


def _encode_local_filename(
    record: Mapping[str, Any],
    accession: str,
    basename: str,
    kind: str,
) -> Optional[str]:
    if basename in SEGWAY_AGGREGATE_FILES:
        return basename
    if kind == "caas":
        return "caas.bed.gz"
    if kind == "label_info":
        return "label_info.tab"
    if kind == "encyclopedia":
        return "segway_encyclopedia.bed.gz"

    submitted_name = _encode_celltype_name_from_submitted_path(record)
    if submitted_name is not None:
        return submitted_name

    biosample_name = _encode_biosample_name(record)
    if biosample_name is not None:
        return _safe_segway_bed_filename(biosample_name)

    if basename and not _is_encode_generic_basename(basename):
        return _ensure_bed_gz(basename)

    return "{}.bed.gz".format(accession)


def _encode_download_url(
    record: Mapping[str, Any],
    accession: str,
    page_url: str,
) -> str:
    for field in ("href", "download_path"):
        value = record.get(field)
        if isinstance(value, str) and value:
            return urllib.parse.urljoin(ENCODE_BASE_URL, value)

    basename = _last_path_name(str(record.get("submitted_file_name", "")))
    if basename and "." in basename:
        extension = ".".join(basename.split(".")[1:])
    else:
        file_format = str(record.get("file_format", "bed"))
        extension = "bed.gz" if file_format == "bed" else file_format
    filename = "{}.{}".format(accession, extension)
    return urllib.parse.urljoin(
        page_url,
        "/files/{}/@@download/{}".format(accession, filename),
    )


def _encode_submitted_basename(record: Mapping[str, Any]) -> str:
    for field in ("submitted_file_name", "href", "download_path"):
        value = record.get(field)
        if isinstance(value, str):
            basename = _last_path_name(value)
            if basename:
                return basename
    accession = _encode_file_accession(record)
    if accession is not None:
        file_format = str(record.get("file_format", "bed"))
        extension = "bed.gz" if file_format == "bed" else file_format
        return "{}.{}".format(accession, extension)
    return ""


def _encode_celltype_name_from_submitted_path(
    record: Mapping[str, Any],
) -> Optional[str]:
    value = record.get("submitted_file_name")
    if not isinstance(value, str):
        return None
    path = urllib.parse.urlparse(value).path or value
    parts = [urllib.parse.unquote(part) for part in Path(path).parts]
    for part in reversed(parts[:-1]):
        name = part.strip().strip("/")
        if _is_encode_celltype_path_component(name):
            return _safe_segway_bed_filename(name)
    return None


def _is_encode_celltype_path_component(value: str) -> bool:
    normalized = value.lower()
    if not normalized:
        return False
    if normalized.startswith("call-"):
        return False
    if normalized in {
        ".",
        "..",
        "files",
        "segway",
        "processed_beds",
        "encode-processing",
        "caper_out_v04_05",
    }:
        return False
    if re.fullmatch(r"[0-9a-f]{8}(-[0-9a-f]{4}){3}-[0-9a-f]{12}", normalized):
        return False
    if re.fullmatch(r"ENC[A-Z]{2}[0-9]{3}[A-Z]{3}", value):
        return False
    return any(character.isalpha() for character in value)


def _encode_biosample_name(record: Mapping[str, Any]) -> Optional[str]:
    for field in ("biosample_term_name", "biosample_summary"):
        value = record.get(field)
        if isinstance(value, str) and value.strip():
            return value

    for field in ("biosample_ontology", "dataset"):
        value = record.get(field)
        name = _encode_nested_biosample_name(value)
        if name is not None:
            return name
    return None


def _encode_nested_biosample_name(value: Any) -> Optional[str]:
    if isinstance(value, Mapping):
        for field in ("term_name", "name", "title", "biosample_summary"):
            nested = value.get(field)
            if isinstance(nested, str) and nested.strip():
                return nested
        for field in ("biosample_ontology", "biosample"):
            nested = _encode_nested_biosample_name(value.get(field))
            if nested is not None:
                return nested
    return None


def _is_encode_hg19_record(record: Mapping[str, Any], kind: str) -> bool:
    if kind in {"label_info"}:
        return True
    assembly = record.get("assembly")
    if assembly is None:
        return True
    if isinstance(assembly, list):
        assemblies = [str(value).lower() for value in assembly]
    else:
        assemblies = [str(assembly).lower()]
    return any(value in {"hg19", "grch37"} for value in assemblies)


def _is_encode_segway_file(
    record: Mapping[str, Any],
    name: str,
    kind: str,
) -> bool:
    text = _encode_record_text(record)
    if name in SEGWAY_AGGREGATE_FILES:
        return True
    if kind in {"caas", "label_info", "encyclopedia"}:
        return True
    if "call-segway_annotate" in text and "recolor" not in text:
        return False
    if "segway" in text:
        return True
    return "semi-automated genome annotation" in text


def _encode_record_text(record: Mapping[str, Any]) -> str:
    values = []  # type: List[str]
    for field in (
        "accession",
        "aliases",
        "submitted_file_name",
        "output_type",
        "file_format",
        "file_format_type",
        "annotation_type",
        "title",
        "description",
        "href",
        "download_path",
    ):
        _append_text_values(values, record.get(field))
    return " ".join(values).lower()


def _append_text_values(values: List[str], value: Any) -> None:
    if isinstance(value, str):
        values.append(value)
    elif isinstance(value, list):
        for item in value:
            _append_text_values(values, item)
    elif isinstance(value, Mapping):
        for item in value.values():
            _append_text_values(values, item)


def _dedupe_entries(entries: Iterable[SegwayFile]) -> Tuple[SegwayFile, ...]:
    deduped = []  # type: List[SegwayFile]
    seen = set()
    for entry in entries:
        key = (entry.name, entry.kind)
        if key in seen:
            continue
        deduped.append(entry)
        seen.add(key)
    return tuple(deduped)


def _encode_related_search_urls(
    payload: Any,
    page_url: str,
) -> Tuple[str, ...]:
    publication_path = _encode_publication_path(payload, page_url)
    if publication_path is None:
        return tuple()
    return (
        _encode_search_url("File", "references", publication_path),
        _encode_search_url("Annotation", "references", publication_path),
    )


def _encode_publication_path(payload: Any, page_url: str) -> Optional[str]:
    if isinstance(payload, Mapping):
        value = payload.get("@id")
        if isinstance(value, str) and value.startswith("/publications/"):
            return _ensure_trailing_slash(value)
    parsed = urllib.parse.urlparse(page_url)
    if parsed.path.startswith("/publications/"):
        return _ensure_trailing_slash(parsed.path)
    return None


def _encode_search_url(object_type: str, link_field: str, link_path: str) -> str:
    params = urllib.parse.urlencode(
        {
            "type": object_type,
            link_field: link_path,
            "format": "json",
            "limit": "all",
            "frame": "embedded",
        }
    )
    return "{}/search/?{}".format(ENCODE_BASE_URL, params)


def _encode_json_url(url: str) -> str:
    parsed = urllib.parse.urlparse(url)
    query = urllib.parse.parse_qsl(parsed.query, keep_blank_values=True)
    keys = {key for key, _value in query}
    if "format" not in keys:
        query.append(("format", "json"))
    if "frame" not in keys and parsed.path.startswith("/publications/"):
        query.append(("frame", "embedded"))
    return urllib.parse.urlunparse(
        parsed._replace(query=urllib.parse.urlencode(query))
    )


def _fetch_encode_json(url: str) -> Any:
    request = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "sjcab-peak2anno-db",
        },
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        text = response.read().decode("utf-8")
    if _is_human_verification_page(text):
        raise ValueError("human verification/captcha page returned")
    return json.loads(text)


def _is_encode_url(url: str) -> bool:
    return urllib.parse.urlparse(url).netloc.lower().endswith("encodeproject.org")


def _last_path_name(value: str) -> str:
    path = urllib.parse.urlparse(value).path or value
    return urllib.parse.unquote(Path(path).name)


def _ensure_bed_gz(name: str) -> str:
    if name.endswith(".bed.gz"):
        return name
    if name.endswith(".bed"):
        return "{}.gz".format(name)
    return "{}.bed.gz".format(_strip_known_suffix(name))


def _safe_segway_bed_filename(value: str) -> str:
    stem = re.sub(r"[^A-Za-z0-9-]+", "_", value).strip("_").upper()
    if not stem:
        stem = "SEGWAY"
    return "{}.bed.gz".format(stem)


def _is_encode_generic_basename(name: str) -> bool:
    return name.lower() in ENCODE_GENERIC_SEGWAY_BASENAMES


def _ensure_trailing_slash(value: str) -> str:
    return value if value.endswith("/") else value + "/"


def _select_files(
    entries: Sequence[SegwayFile],
    names: Optional[Union[str, Iterable[str]]],
    tissue: Optional[Union[str, Iterable[str]]],
    cellline: Optional[Union[str, Iterable[str]]],
    all_celltypes: bool,
    include_encyclopedia: bool,
    include_caas: bool,
    include_label_info: bool,
) -> Tuple[SegwayFile, ...]:
    explicit_names = _normalize_names(names)
    tissue_queries = _normalize_queries(tissue)
    cellline_queries = _normalize_queries(cellline)
    selected = []  # type: List[SegwayFile]
    seen = set()

    by_alias = _entry_aliases(entries)
    for name in explicit_names:
        key = _normalize_text(_strip_known_suffix(name))
        entry = by_alias.get(key)
        if entry is None:
            raise UnknownResourceError("No Segway file matched {!r}.".format(name))
        _append_selected(selected, seen, entry)

    for filename, include in (
        ("segway_encyclopedia.bed.gz", include_encyclopedia),
        ("caas.bed.gz", include_caas),
        ("label_info.tab", include_label_info),
    ):
        if include:
            entry = by_alias.get(_normalize_text(_strip_known_suffix(filename)))
            if entry is not None:
                _append_selected(selected, seen, entry)

    if all_celltypes:
        for entry in entries:
            if entry.kind == "celltype":
                _append_selected(selected, seen, entry)

    if tissue_queries or cellline_queries:
        queries = tissue_queries + cellline_queries
        for entry in entries:
            if entry.kind == "celltype" and _matches_any(entry.search_text, queries):
                _append_selected(selected, seen, entry)

    if not selected and not (
        explicit_names
        or tissue_queries
        or cellline_queries
        or all_celltypes
        or include_encyclopedia
        or include_caas
        or include_label_info
    ):
        entry = by_alias.get(_normalize_text(_strip_known_suffix(SEGWAY_DEFAULT_FILE)))
        if entry is not None:
            _append_selected(selected, seen, entry)

    return tuple(selected)


def _write_segway_metadata_tsv(
    outputs: Mapping[str, Path],
    entries: Mapping[str, SegwayFile],
    genome: str,
    data_dir: Optional[PathLike],
    progress: Optional[ProgressCallback],
) -> Path:
    report_progress(progress, "download-segway: writing metadata TSV")
    output = segway_metadata_tsv_path(data_dir, genome)
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output.with_name(output.name + ".tmp")
    rows = _read_existing_metadata_rows(output)
    for key, destination in outputs.items():
        entry = entries[key]
        rows[_metadata_row_key(destination, genome, entry.key)] = (
            str(destination.expanduser()),
            genome,
            entry.key,
            entry.kind,
            entry.encode_url,
            entry.url,
        )
    try:
        with tmp_path.open("w", encoding="utf-8") as handle:
            handle.write("{}\n".format("\t".join(SEGWAY_METADATA_FIELDS)))
            for key in sorted(rows):
                handle.write("{}\n".format("\t".join(rows[key])))
        tmp_path.replace(output)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()
    report_progress(progress, "download-segway: metadata TSV done")
    return output


def _install_segway_download_manifest(data_dir: Optional[PathLike]) -> Path:
    source = _bundled_segway_manifest_path()
    output = segway_download_manifest_path(data_dir)
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output.with_name(output.name + ".tmp")
    try:
        tmp_path.write_bytes(source.read_bytes())
        tmp_path.replace(output)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()
    return output


def _bundled_segway_manifest_path() -> Path:
    return Path(__file__).with_name("data") / SEGWAY_DOWNLOAD_MANIFEST_NAME


def _read_existing_metadata_rows(path: Path) -> Dict[Tuple[str, str, str], Tuple[str, ...]]:
    if not path.exists():
        return {}

    rows = {}  # type: Dict[Tuple[str, str, str], Tuple[str, ...]]
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if tuple(header) != tuple(SEGWAY_METADATA_FIELDS):
            return {}
        for line in handle:
            values = tuple(line.rstrip("\n").split("\t"))
            if len(values) != len(SEGWAY_METADATA_FIELDS):
                continue
            rows[(values[0], values[1], values[2])] = values
    return rows


def _metadata_row_key(destination: Path, genome: str, name: str) -> Tuple[str, str, str]:
    return (str(destination.expanduser()), genome, name)


def _append_entry(
    entries: List[SegwayFile],
    seen: set,
    name: str,
    url: str,
    kind: str,
    encode_url: str = "",
) -> None:
    key = (name, kind)
    if key in seen:
        return
    entries.append(SegwayFile(name=name, url=url, kind=kind, encode_url=encode_url))
    seen.add(key)


def _append_selected(
    selected: List[SegwayFile],
    seen: set,
    entry: SegwayFile,
) -> None:
    if entry.name in seen:
        return
    selected.append(entry)
    seen.add(entry.name)


def _entry_aliases(entries: Sequence[SegwayFile]) -> Dict[str, SegwayFile]:
    aliases = {}  # type: Dict[str, SegwayFile]
    for entry in entries:
        values = [
            entry.name,
            _strip_known_suffix(entry.name),
            entry.key,
        ]
        if entry.name == "segway_encyclopedia.bed.gz":
            values.extend(["encyclopedia", "segway encyclopedia"])
        elif entry.name == "caas.bed.gz":
            values.extend(["caas", "position-wise caas", "positionwise caas"])
        elif entry.name == "label_info.tab":
            values.extend(["label info", "labels", "label_info"])
        for value in values:
            aliases[_normalize_text(value)] = entry
    return aliases


def _interpreted_directory_urls(html: str, page_url: str) -> Tuple[str, ...]:
    urls = []
    seen = set()
    for href in _html_links(html):
        normalized = href.rstrip("/")
        if normalized.endswith("interpreted"):
            url = urllib.parse.urljoin(page_url, normalized + "/")
            if url not in seen:
                urls.append(url)
                seen.add(url)
    default_url = urllib.parse.urljoin(page_url, SEGWAY_INTERPRETED_DIR)
    if default_url not in seen:
        urls.append(default_url)
    return tuple(urls)


def _html_links(html: str) -> Tuple[str, ...]:
    return tuple(
        match.group(1)
        for match in re.finditer(
            r"""href\s*=\s*["']?([^"'\s>]+)""",
            html,
            re.IGNORECASE,
        )
    )


def _is_human_verification_page(html: str) -> bool:
    normalized = html.lower()
    return (
        "human verification" in normalized
        or "captcha" in normalized
        or "x-amzn-waf-action" in normalized
    )


def _aggregate_kind(filename: str) -> str:
    if filename == "segway_encyclopedia.bed.gz":
        return "encyclopedia"
    if filename == "caas.bed.gz":
        return "caas"
    if filename == "label_info.tab":
        return "label_info"
    return "aggregate"


def _normalize_genome(genome: str) -> str:
    value = str(genome).lower()
    if value not in SEGWAY_GENOMES:
        raise UnknownResourceError(
            "Segway encyclopedia downloads are hg19 only. For {}, first download "
            "hg19 and run a liftover script generated by install-segway or "
            "download-segway with --yes-liftover.".format(genome)
        )
    return value


def _normalize_target_genome(genome: str) -> str:
    value = str(genome).strip()
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]*", value):
        raise UnknownResourceError(
            "Target genome must be a simple UCSC-style genome build such as hg38."
        )
    if value.lower() == "hg19":
        raise UnknownResourceError("Segway resources are already hg19.")
    return value


def _ucsc_target_token(genome: str) -> str:
    return genome[:1].upper() + genome[1:]


def _normalize_names(
    names: Optional[Union[str, Iterable[str]]]
) -> Tuple[str, ...]:
    if names is None:
        return tuple()
    raw_values = []  # type: List[str]
    if isinstance(names, str):
        raw_values.extend(names.split(","))
    else:
        for value in names:
            raw_values.extend(str(value).split(","))
    return tuple(value.strip() for value in raw_values if value.strip())


def _normalize_queries(
    queries: Optional[Union[str, Iterable[str]]]
) -> Tuple[str, ...]:
    if queries is None:
        return tuple()
    if isinstance(queries, str):
        values = [queries]
    else:
        values = [str(value) for value in queries]
    return tuple(value for value in values if value.strip())


def _matches_any(text: str, queries: Tuple[str, ...]) -> bool:
    if not queries:
        return False
    normalized_text = _normalize_text(text)
    for query in queries:
        normalized_query = _normalize_text(query)
        query_tokens = [token for token in normalized_query.split() if token]
        if query_tokens and all(token in normalized_text for token in query_tokens):
            return True
        if SequenceMatcher(None, normalized_query, normalized_text).ratio() >= 0.72:
            return True
    return False


def _normalize_text(value: str) -> str:
    return " ".join(re.sub(r"[^a-z0-9]+", " ", value.lower()).split())


def _strip_known_suffix(name: str) -> str:
    value = name.strip()
    for suffix in (".bed.gz", ".tab"):
        if value.endswith(suffix):
            return value[: -len(suffix)]
    return value


def _shell_quote(value: str) -> str:
    return "'{}'".format(value.replace("'", "'\"'\"'"))


def _fetch_text(url: str) -> str:
    with urllib.request.urlopen(url, timeout=60) as response:
        return response.read().decode("utf-8")


def _download_file(
    url: str,
    destination: Path,
    progress: Optional[ProgressCallback] = None,
) -> None:
    download_file(url, destination, timeout=120, progress=progress)
