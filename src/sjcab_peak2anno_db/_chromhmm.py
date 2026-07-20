"""Download Roadmap Epigenomics ChromHMM dense BED resources."""

from __future__ import annotations

import json
import re
import urllib.request
from dataclasses import dataclass
from difflib import SequenceMatcher
from pathlib import Path
from typing import Dict, Iterable, Mapping, Optional, Tuple, Union

from ._download import ProgressCallback, download_file, report_progress
from ._download_log import record_download_url
from ._registry import UnknownResourceError, user_data_dir

PathLike = Union[str, Path]

CHROMHMM_DIR_NAME = "chromhmm"
CHROMHMM_GENOMES = ("hg19", "hg38")
CHROMHMM_METADATA_URL = "https://egg2.wustl.edu/roadmap/web_portal/data.js"
CHROMHMM_METADATA_TSV_NAME = "metadata.tsv"
CHROMHMM_METADATA_FIELDS = (
    "path",
    "version",
    "state_number",
    "eid",
    "group",
    "mnemonic",
    "name",
    "name_edacc_r9",
    "category",
)
_MODEL_CONFIGS = {
    15: (
        "coreMarks",
        "https://egg2.wustl.edu/roadmap/data/byFileType/"
        "chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/",
    ),
    18: (
        "core_K27ac",
        "https://egg2.wustl.edu/roadmap/data/byFileType/"
        "chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/",
    ),
    25: (
        "imputed12marks",
        "https://egg2.wustl.edu/roadmap/data/byFileType/"
        "chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/",
    ),
}


@dataclass(frozen=True)
class ChromHMMEpigenome:
    """Roadmap epigenome metadata used for ChromHMM sample selection."""

    eid: str
    group: str
    mnemonic: str
    name: str
    name_edacc_r9: str
    category: str

    @property
    def tissue_text(self) -> str:
        return " ".join((self.group, self.mnemonic, self.name, self.category))

    @property
    def cellline_text(self) -> str:
        return " ".join((self.name, self.name_edacc_r9, self.mnemonic))


def chromhmm_dir(
    data_dir: Optional[PathLike] = None,
    model: int = 18,
    genome: str = "hg19",
) -> Path:
    """Return the local cache directory for one ChromHMM model."""

    model = _normalize_model(model)
    genome = _normalize_genome(genome)
    return chromhmm_root(data_dir) / genome / "{}state".format(model)


def chromhmm_root(data_dir: Optional[PathLike] = None) -> Path:
    """Return the local ChromHMM cache root."""

    return user_data_dir(data_dir) / CHROMHMM_DIR_NAME


def chromhmm_filename(eid: str, model: int = 18) -> str:
    """Return the Roadmap dense BED filename for one epigenome/model pair."""

    eid = _normalize_id(eid)
    model = _normalize_model(model)
    model_label = _MODEL_CONFIGS[model][0]
    return "{}_{}_{}_dense.bed.gz".format(eid, model, model_label)


def chromhmm_path(
    eid: str,
    data_dir: Optional[PathLike] = None,
    model: int = 18,
    genome: str = "hg19",
) -> Path:
    """Return the local cache path for one ChromHMM dense BED file."""

    return chromhmm_dir(data_dir, model, genome) / chromhmm_filename(eid, model)


def chromhmm_metadata_tsv_path(
    data_dir: Optional[PathLike] = None,
    model: int = 18,
    genome: str = "hg19",
) -> Path:
    """Return the local global ChromHMM metadata TSV path."""

    _normalize_model(model)
    _normalize_genome(genome)
    return chromhmm_root(data_dir) / CHROMHMM_METADATA_TSV_NAME


def chromhmm_url(eid: str, model: int = 18, genome: str = "hg19") -> str:
    """Return the Roadmap dense BED URL for one epigenome/model pair."""

    model = _normalize_model(model)
    return _chromhmm_model_url(model, genome) + chromhmm_filename(eid, model)


def chromhmm_metadata(
    metadata_url: str = CHROMHMM_METADATA_URL,
) -> Tuple[ChromHMMEpigenome, ...]:
    """Return Roadmap epigenome metadata from the portal's ``data.js``."""

    return _parse_metadata(_fetch_text(metadata_url))


def chromhmm_available_ids(model: int = 18, genome: str = "hg19") -> Tuple[str, ...]:
    """Return epigenome IDs with dense BED files listed for one model."""

    model = _normalize_model(model)
    html = _fetch_text(_chromhmm_model_url(model, genome))
    pattern = re.compile(
        r"\b(E\d{{3}})_{}_{}_dense\.bed\.gz\b".format(
            model, re.escape(_MODEL_CONFIGS[model][0])
        )
    )
    ids = []
    seen = set()
    for match in pattern.finditer(html):
        eid = match.group(1)
        if eid in seen:
            continue
        seen.add(eid)
        ids.append(eid)
    return tuple(ids)


def download_chromhmm(
    data_dir: Optional[PathLike] = None,
    model: int = 18,
    genome: str = "hg19",
    ids: Optional[Union[str, Iterable[str]]] = None,
    tissue: Optional[Union[str, Iterable[str]]] = None,
    cellline: Optional[Union[str, Iterable[str]]] = None,
    overwrite: bool = True,
    metadata_url: str = CHROMHMM_METADATA_URL,
    progress: Optional[ProgressCallback] = None,
) -> Mapping[str, Path]:
    """Download Roadmap ChromHMM dense BED files.

    ``model`` can be ``15``, ``18`` (default), or ``25``. ``genome`` can be
    ``hg19`` (default) or ``hg38`` for lifted-over Roadmap BED files. Select
    epigenomes by explicit IDs such as ``"E001,E063"``, tissue/group keywords,
    or cell-line keywords. With no selector, all dense BED files available for
    the selected model/genome are downloaded.
    """

    model = _normalize_model(model)
    genome = _normalize_genome(genome)
    selected = _select_ids(
        model,
        genome,
        ids=ids,
        tissue=tissue,
        cellline=cellline,
        metadata_url=metadata_url,
    )
    if not selected:
        raise UnknownResourceError("No ChromHMM epigenomes matched the request.")

    target_dir = chromhmm_dir(data_dir, model, genome)
    target_dir.mkdir(parents=True, exist_ok=True)

    outputs = {}  # type: Dict[str, Path]
    for eid in selected:
        destination = target_dir / chromhmm_filename(eid, model)
        url = chromhmm_url(eid, model, genome)
        if overwrite or not destination.exists():
            if progress is None:
                _download_file(url, destination)
            else:
                _download_file(url, destination, progress=progress)
            record_download_url(url, data_dir=data_dir, destination=destination)
        outputs[eid] = destination
    _write_chromhmm_metadata_tsv(
        outputs,
        genome=genome,
        model=model,
        data_dir=data_dir,
        metadata_url=metadata_url,
        progress=progress,
    )
    return outputs


def _select_ids(
    model: int,
    genome: str,
    ids: Optional[Union[str, Iterable[str]]],
    tissue: Optional[Union[str, Iterable[str]]],
    cellline: Optional[Union[str, Iterable[str]]],
    metadata_url: str,
) -> Tuple[str, ...]:
    explicit_ids = list(_normalize_ids(ids))
    tissue_queries = _normalize_queries(tissue)
    cellline_queries = _normalize_queries(cellline)

    available = tuple()
    try:
        available = chromhmm_available_ids(model, genome)
    except Exception:
        available = tuple()
    available_set = set(available)

    selected = []
    seen = set()
    for eid in explicit_ids:
        if available_set and eid not in available_set:
            continue
        selected.append(eid)
        seen.add(eid)

    if tissue_queries or cellline_queries:
        metadata = chromhmm_metadata(metadata_url)
        for entry in metadata:
            if available_set and entry.eid not in available_set:
                continue
            if (
                _matches_any(entry.tissue_text, tissue_queries)
                or _matches_any(entry.cellline_text, cellline_queries)
            ) and entry.eid not in seen:
                selected.append(entry.eid)
                seen.add(entry.eid)

    if not selected and not explicit_ids and not tissue_queries and not cellline_queries:
        if available:
            selected = list(available)
        else:
            selected = [entry.eid for entry in chromhmm_metadata(metadata_url)]

    return tuple(selected)


def _parse_metadata(text: str) -> Tuple[ChromHMMEpigenome, ...]:
    match = re.search(r"var\s+data_epg\s*=\s*(\[.*?\]);", text, re.DOTALL)
    if not match:
        raise ValueError("Could not find Roadmap data_epg metadata in data.js.")

    payload = re.sub(r",\s*([}\]])", r"\1", match.group(1))
    values = json.loads(payload)
    entries = []
    for value in values:
        eid = value.get("eid", "")
        if not eid:
            continue
        entries.append(
            ChromHMMEpigenome(
                eid=eid,
                group=value.get("group", ""),
                mnemonic=value.get("mnemonic", ""),
                name=value.get("name", ""),
                name_edacc_r9=value.get("name_EDACC_R9", ""),
                category=value.get("class", ""),
            )
        )
    return tuple(entries)


def _write_chromhmm_metadata_tsv(
    outputs: Mapping[str, Path],
    genome: str,
    model: int,
    data_dir: Optional[PathLike],
    metadata_url: str,
    progress: Optional[ProgressCallback],
) -> Path:
    report_progress(progress, "download-chromhmm: writing metadata TSV")
    metadata_by_id: Dict[str, ChromHMMEpigenome] = {}
    try:
        metadata_by_id = {entry.eid: entry for entry in chromhmm_metadata(metadata_url)}
    except Exception:
        metadata_by_id = {}

    output = chromhmm_metadata_tsv_path(data_dir)
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output.with_name(output.name + ".tmp")
    rows = _read_existing_metadata_rows(output)
    for eid, destination in outputs.items():
        rows[_metadata_row_key(destination, genome, model)] = _metadata_row(
            eid,
            destination,
            genome,
            model,
            metadata_by_id.get(eid),
        )
    try:
        with tmp_path.open("w", encoding="utf-8") as handle:
            handle.write("{}\n".format("\t".join(CHROMHMM_METADATA_FIELDS)))
            for key in sorted(rows):
                handle.write("{}\n".format("\t".join(rows[key])))
        tmp_path.replace(output)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()
    if metadata_by_id:
        record_download_url(metadata_url, data_dir=data_dir, destination=output)
    report_progress(progress, "download-chromhmm: metadata TSV done")
    return output


def _read_existing_metadata_rows(path: Path) -> Dict[Tuple[str, str, str], Tuple[str, ...]]:
    if not path.exists():
        return {}

    rows = {}  # type: Dict[Tuple[str, str, str], Tuple[str, ...]]
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if tuple(header) != CHROMHMM_METADATA_FIELDS:
            return {}
        for line in handle:
            values = tuple(line.rstrip("\n").split("\t"))
            if len(values) != len(CHROMHMM_METADATA_FIELDS):
                continue
            rows[(values[0], values[1], values[2])] = values
    return rows


def _metadata_row_key(
    destination: Path,
    genome: str,
    model: int,
) -> Tuple[str, str, str]:
    return (str(destination.expanduser()), genome, str(model))


def _metadata_row(
    eid: str,
    destination: Path,
    genome: str,
    model: int,
    entry: Optional[ChromHMMEpigenome],
) -> Tuple[str, ...]:
    if entry is None:
        metadata_values = (eid, "", "", "", "", "")
    else:
        metadata_values = (
            entry.eid,
            entry.group,
            entry.mnemonic,
            entry.name,
            entry.name_edacc_r9,
            entry.category,
        )
    return (
        str(destination.expanduser()),
        genome,
        str(model),
    ) + tuple(_sanitize_metadata_value(value) for value in metadata_values)


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


def _normalize_model(model: int) -> int:
    try:
        value = int(model)
    except (TypeError, ValueError) as exc:
        raise UnknownResourceError("ChromHMM model must be 15, 18, or 25.") from exc
    if value not in _MODEL_CONFIGS:
        raise UnknownResourceError("ChromHMM model must be 15, 18, or 25.")
    return value


def _normalize_genome(genome: str) -> str:
    value = str(genome).lower()
    if value not in CHROMHMM_GENOMES:
        raise UnknownResourceError(
            "ChromHMM genome must be hg19 or hg38."
        )
    return value


def _chromhmm_model_url(model: int, genome: str = "hg19") -> str:
    model = _normalize_model(model)
    genome = _normalize_genome(genome)
    base = _MODEL_CONFIGS[model][1]
    if genome == "hg38":
        return base + "bed_hg38_lifted_over/"
    return base


def _normalize_ids(ids: Optional[Union[str, Iterable[str]]]) -> Tuple[str, ...]:
    if ids is None:
        return tuple()
    if isinstance(ids, str):
        raw_values = re.split(r"[\s,]+", ids)
    else:
        raw_values = []
        for value in ids:
            raw_values.extend(re.split(r"[\s,]+", str(value)))

    normalized = []
    for value in raw_values:
        if not value:
            continue
        normalized.append(_normalize_id(value))
    return tuple(normalized)


def _normalize_id(value: str) -> str:
    eid = value.strip().upper()
    if not re.fullmatch(r"E\d{3}", eid):
        raise UnknownResourceError(
            "ChromHMM epigenome IDs should look like E001 or E063."
        )
    return eid


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


def _normalize_text(value: str) -> str:
    return " ".join(re.sub(r"[^a-z0-9]+", " ", value.lower()).split())


def _sanitize_metadata_value(value: str) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")
    return sanitized or "."


def _fetch_text(url: str) -> str:
    with urllib.request.urlopen(url, timeout=60) as response:
        return response.read().decode("utf-8")


def _download_file(
    url: str,
    destination: Path,
    progress: Optional[ProgressCallback] = None,
) -> None:
    download_file(url, destination, timeout=120, progress=progress)
