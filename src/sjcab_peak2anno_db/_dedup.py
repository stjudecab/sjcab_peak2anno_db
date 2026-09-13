"""Select or filter one GENCODE isoform per gene."""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Mapping, Optional, Sequence, Tuple, Union

from ._derive import _length, _parse_bed_fields, _site_interval, write_tes, write_tss
from ._config import (
    configured_bed_score_column,
    configured_peak_txt_delimiter,
    configured_txt_score_column,
)
from ._registry import default_version as registry_default_version
from ._registry import path as registry_path

PathLike = Union[str, os.PathLike]

DEDUP_METHODS = ("longcol5", "long", "peak", "isoID", "isoexp", "perover")
@dataclass(frozen=True)
class _IsoformRecord:
    fields: Tuple[str, ...]
    line: str
    chrom: str
    start: int
    end: int
    gene_key: str
    transcript_id: str
    transcript_id_base: str
    strand: str
    length: float
    order: int


def dedup_gencode_bed(
    method: str = "longcol5",
    selector: Optional[PathLike] = None,
    output_dir: PathLike = ".",
    gene_bed: Optional[PathLike] = None,
    species: Optional[str] = None,
    version: Optional[str] = "default",
    data_dir: Optional[PathLike] = None,
    promoter_bp: Union[int, str] = 2000,
    inclusive: bool = True,
    output_prefix: Optional[str] = None,
    gene_key: str = "symbol",
    promoter_down_bp: Optional[Union[int, str]] = None,
) -> Mapping[str, Path]:
    """Keep one isoform per gene based on a selector file.

    ``method`` is one of ``longcol5``, ``long``, ``peak``, ``isoID``,
    ``isoexp``, or ``perover``. ``longcol5`` is the default length selector;
    it uses BED column 5. ``long`` uses ``end - start`` and needs no selector.
    Genes with no selected/overlapping isoform fall back to the longest isoform.
    Use :func:`filter_gencode_bed` to omit those genes instead.
    """

    return _select_gencode_bed(
        method=method,
        selector=selector,
        output_dir=output_dir,
        gene_bed=gene_bed,
        species=species,
        version=version,
        data_dir=data_dir,
        promoter_bp=promoter_bp,
        promoter_down_bp=promoter_down_bp,
        inclusive=inclusive,
        output_prefix=output_prefix,
        gene_key=gene_key,
        fallback_to_longest=True,
        default_prefix_kind="dedup",
    )


def filter_gencode_bed(
    method: str,
    selector: Optional[PathLike] = None,
    output_dir: PathLike = ".",
    gene_bed: Optional[PathLike] = None,
    species: Optional[str] = None,
    version: Optional[str] = "default",
    data_dir: Optional[PathLike] = None,
    promoter_bp: Union[int, str] = 2000,
    inclusive: bool = True,
    output_prefix: Optional[str] = None,
    gene_key: str = "symbol",
    promoter_down_bp: Optional[Union[int, str]] = None,
) -> Mapping[str, Path]:
    """Keep one isoform per gene and omit genes with no selector support."""

    return _select_gencode_bed(
        method=method,
        selector=selector,
        output_dir=output_dir,
        gene_bed=gene_bed,
        species=species,
        version=version,
        data_dir=data_dir,
        promoter_bp=promoter_bp,
        promoter_down_bp=promoter_down_bp,
        inclusive=inclusive,
        output_prefix=output_prefix,
        gene_key=gene_key,
        fallback_to_longest=False,
        default_prefix_kind="filter",
    )


def _select_gencode_bed(
    method: str,
    selector: PathLike,
    output_dir: PathLike,
    gene_bed: Optional[PathLike],
    species: Optional[str],
    version: Optional[str],
    data_dir: Optional[PathLike],
    promoter_bp: Union[int, str],
    promoter_down_bp: Optional[Union[int, str]],
    inclusive: bool,
    output_prefix: Optional[str],
    gene_key: str,
    fallback_to_longest: bool,
    default_prefix_kind: str,
) -> Mapping[str, Path]:
    selected_method = _normalize_method(method)
    key_mode = _normalize_gene_key(gene_key)
    gene_bed_path, resolved_species, resolved_version = _resolve_gene_bed(
        gene_bed,
        species,
        version,
        data_dir,
    )
    records = _read_gene_bed(gene_bed_path, key_mode)
    records_by_gene = _group_by_gene(records)
    promoter = _parse_bp(promoter_bp)
    promoter_down = _parse_bp(
        promoter_bp if promoter_down_bp is None else promoter_down_bp
    )

    if selected_method == "longcol5":
        scores = _score_length(records, use_column_five=True)
    elif selected_method == "long":
        scores = _score_length(records, use_column_five=False)
    elif selected_method == "peak":
        if selector is None:
            raise ValueError("peak requires a selector file.")
        scores = _score_peak(records, selector, promoter, promoter_down, inclusive)
    elif selected_method == "perover":
        if selector is None:
            raise ValueError("perover requires a selector file.")
        scores = _score_perover(records, selector, promoter, promoter_down, inclusive)
    elif selected_method == "isoID":
        if selector is None:
            raise ValueError("isoID requires a selector file.")
        scores = _score_iso_ids(records, selector, inclusive)
    else:
        if selector is None:
            raise ValueError("isoexp requires a selector file.")
        scores = _score_iso_expression(records, selector, inclusive)

    selected = []
    for gene_key in sorted(records_by_gene, key=lambda key: records_by_gene[key][0].order):
        record = _select_record(records_by_gene[gene_key], scores, fallback_to_longest)
        if record is not None:
            selected.append(record)

    prefix = output_prefix or _default_output_prefix(
        gene_bed_path,
        resolved_species,
        resolved_version,
        default_prefix_kind,
        selected_method,
    )
    output = Path(output_dir).expanduser()
    output.mkdir(parents=True, exist_ok=True)
    gene_output = output / "{}.gene.bed".format(prefix)
    with gene_output.open("w", encoding="utf-8") as handle:
        for record in selected:
            handle.write(record.line)
            if not record.line.endswith("\n"):
                handle.write("\n")

    tss_output = output / "{}.tss.bed".format(prefix)
    tes_output = output / "{}.tes.bed".format(prefix)
    write_tss(gene_output, tss_output)
    write_tes(gene_output, tes_output)
    return {"gene": gene_output, "tss": tss_output, "tes": tes_output}


def _resolve_gene_bed(
    gene_bed: Optional[PathLike],
    species: Optional[str],
    version: Optional[str],
    data_dir: Optional[PathLike],
) -> Tuple[Path, Optional[str], Optional[str]]:
    if gene_bed is not None:
        return Path(gene_bed).expanduser(), None, None
    if not species:
        raise ValueError("Provide gene_bed or species/version.")
    selected_version = version or "default"
    resolved_version = (
        registry_default_version(species, "gene", "all")
        if selected_version in {"default", "def", "latest"}
        else selected_version
    )
    return registry_path(
        species,
        "gene",
        selected_version,
        data_dir=data_dir,
        prefer_user=True,
        isoform_set="all",
    ), species, resolved_version


def _read_gene_bed(
    gene_bed: Path,
    gene_key: str,
) -> List[_IsoformRecord]:
    records = []
    with gene_bed.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = tuple(_parse_bed_fields(line, line_number))
            transcript_id = _field_or_empty(fields, 7)
            record_gene_key = _record_gene_key(fields, gene_key)
            records.append(
                _IsoformRecord(
                    fields=fields,
                    line=line,
                    chrom=fields[0],
                    start=int(fields[1]),
                    end=int(fields[2]),
                    gene_key=record_gene_key,
                    transcript_id=transcript_id,
                    transcript_id_base=_field_or_empty(fields, 10)
                    or _strip_version(transcript_id),
                    strand=fields[5],
                    length=_length(list(fields)),
                    order=len(records),
                )
            )
    return records


def _group_by_gene(
    records: Sequence[_IsoformRecord],
) -> Mapping[str, List[_IsoformRecord]]:
    grouped = {}  # type: Dict[str, List[_IsoformRecord]]
    for record in records:
        grouped.setdefault(record.gene_key, []).append(record)
    return grouped


def _record_gene_key(
    fields: Sequence[str],
    gene_key: str,
) -> str:
    if gene_key == "ensid":
        gene_id = _field_or_empty(fields, 6)
        if gene_id:
            return _strip_version(gene_id)
    return _field_or_empty(fields, 3)


def _default_output_prefix(
    gene_bed: Path,
    species: Optional[str],
    version: Optional[str],
    prefix_kind: str,
    method: str,
) -> str:
    if species and version:
        base = "{}.{}".format(species, version)
    else:
        base = _input_bed_prefix(gene_bed)
    return "{}.{}{}".format(base, prefix_kind, method)


def _input_bed_prefix(gene_bed: Path) -> str:
    name = gene_bed.name
    if name.endswith(".bed"):
        return name[:-4]
    return gene_bed.stem


def _select_record(
    records: Sequence[_IsoformRecord],
    scores: Mapping[int, float],
    fallback_to_longest: bool,
) -> Optional[_IsoformRecord]:
    scored = [record for record in records if record.order in scores]
    if scored:
        return max(
            scored,
            key=lambda record: (scores[record.order], record.length, -record.order),
        )
    if not fallback_to_longest:
        return None
    return max(records, key=lambda record: (record.length, -record.order))


def _score_iso_ids(
    records: Sequence[_IsoformRecord],
    selector: PathLike,
    inclusive: bool,
) -> Mapping[int, float]:
    ranked_ids = _read_ranked_ids(selector)
    scores = {}  # type: Dict[int, float]
    for record in records:
        matched_rank = _match_identifier(record, ranked_ids, inclusive)
        if matched_rank is not None:
            scores[record.order] = -float(matched_rank)
    return scores


def _score_length(
    records: Sequence[_IsoformRecord],
    use_column_five: bool,
) -> Mapping[int, float]:
    scores = {}  # type: Dict[int, float]
    for record in records:
        if use_column_five:
            try:
                scores[record.order] = float(record.fields[4])
            except (IndexError, ValueError) as exc:
                raise ValueError(
                    "longcol5 requires a numeric BED column 5 for every isoform."
                ) from exc
        else:
            scores[record.order] = float(record.end - record.start)
    return scores


def _score_iso_expression(
    records: Sequence[_IsoformRecord],
    selector: PathLike,
    inclusive: bool,
) -> Mapping[int, float]:
    expression = _read_expression_table(selector)
    scores = {}  # type: Dict[int, float]
    for record in records:
        values = []
        for identifier in _record_identifiers(record, inclusive):
            if identifier in expression:
                values.append(expression[identifier])
        if values:
            scores[record.order] = max(values)
    return scores


def _score_peak(
    records: Sequence[_IsoformRecord],
    selector: PathLike,
    promoter_bp: int,
    promoter_down_bp: int,
    inclusive: bool,
) -> Mapping[int, float]:
    intervals_by_chrom = _read_scored_bed(selector)
    scores = {}  # type: Dict[int, float]
    for record in records:
        promoter = _promoter_interval(record, promoter_bp, promoter_down_bp)
        values = [
            score
            for start, end, score in _candidate_intervals(
                intervals_by_chrom, record.chrom, promoter
            )
            if _overlaps(promoter, (start, end), inclusive)
        ]
        if values:
            scores[record.order] = max(values)
    return scores


def _score_perover(
    records: Sequence[_IsoformRecord],
    selector: PathLike,
    promoter_bp: int,
    promoter_down_bp: int,
    inclusive: bool,
) -> Mapping[int, float]:
    intervals_by_chrom = _read_bed(selector)
    scores = {}  # type: Dict[int, float]
    for record in records:
        promoter = _promoter_interval(record, promoter_bp, promoter_down_bp)
        pieces = []
        for start, end in _candidate_intervals(intervals_by_chrom, record.chrom, promoter):
            if not _overlaps(promoter, (start, end), inclusive):
                continue
            pieces.append((max(promoter[0], start), min(promoter[1], end)))
        if not pieces:
            continue
        overlap_bp = sum(end - start for start, end in _merge_intervals(pieces))
        promoter_bp_value = max(1, promoter[1] - promoter[0])
        scores[record.order] = overlap_bp / float(promoter_bp_value)
    return scores


def _candidate_intervals(
    intervals_by_chrom: Mapping[str, Sequence[Tuple]],
    chrom: str,
    query: Tuple[int, int],
) -> Iterable[Tuple]:
    query_start, query_end = query
    for interval in intervals_by_chrom.get(chrom, ()):
        start, end = interval[0], interval[1]
        if end <= query_start:
            continue
        if start >= query_end:
            break
        yield interval


def _overlaps(
    query: Tuple[int, int],
    interval: Tuple[int, int],
    inclusive: bool,
) -> bool:
    query_start, query_end = query
    interval_start, interval_end = interval
    if inclusive:
        return max(query_start, interval_start) < min(query_end, interval_end)
    return interval_start >= query_start and interval_end <= query_end


def _promoter_interval(
    record: _IsoformRecord,
    promoter_bp: int,
    promoter_down_bp: int,
) -> Tuple[int, int]:
    tss = int(_site_interval(list(record.fields), "tss")[0])
    if record.strand == "-":
        return max(0, tss - promoter_down_bp), tss + promoter_bp + 1
    return max(0, tss - promoter_bp), tss + promoter_down_bp + 1


def _read_scored_bed(
    bed_path: PathLike,
) -> Mapping[str, List[Tuple[int, int, float]]]:
    scored = {}  # type: Dict[str, List[Tuple[int, int, float]]]
    for line_number, fields in _iter_selector_fields(bed_path):
        try:
            chrom, start, end = _selector_coordinates(fields, bed_path, line_number)
            if _is_bed_row(fields):
                score = _column_score(
                    fields, configured_bed_score_column(), bed_path, line_number
                )
            else:
                score = _column_score(
                    fields, configured_txt_score_column(), bed_path, line_number
                )
        except (IndexError, ValueError) as exc:
            raise ValueError(
                "{}:{} has invalid interval or peak score.".format(
                    bed_path, line_number
                )
            ) from exc
        if end > start:
            scored.setdefault(chrom, []).append((start, end, score))
    return _sort_interval_map(scored)


def _read_bed(
    bed_path: PathLike,
) -> Mapping[str, List[Tuple[int, int]]]:
    intervals = {}  # type: Dict[str, List[Tuple[int, int]]]
    for line_number, fields in _iter_selector_fields(bed_path):
        try:
            chrom, start, end = _selector_coordinates(fields, bed_path, line_number)
        except (IndexError, ValueError) as exc:
            raise ValueError(
                "{}:{} has invalid interval.".format(bed_path, line_number)
            ) from exc
        if end > start:
            intervals.setdefault(chrom, []).append((start, end))
    return _sort_interval_map(intervals)


def _iter_selector_fields(selector: PathLike) -> Iterable[Tuple[int, List[str]]]:
    with Path(selector).expanduser().open("r", encoding="utf-8") as handle:
        header_skipped = False
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = _split_selector_fields(line)
            if not header_skipped and not _looks_like_interval(fields):
                header_skipped = True
                continue
            header_skipped = True
            if _looks_like_interval(fields):
                yield line_number, fields


def _looks_like_interval(fields: Sequence[str]) -> bool:
    if len(fields) >= 3:
        try:
            int(fields[1])
            int(fields[2])
            return True
        except ValueError:
            pass
    return bool(fields) and _region_match(fields[0]) is not None


def _is_bed_row(fields: Sequence[str]) -> bool:
    if len(fields) < 3:
        return False
    try:
        int(fields[1])
        int(fields[2])
    except ValueError:
        return False
    return True


def _selector_coordinates(
    fields: Sequence[str],
    selector: PathLike,
    line_number: int,
) -> Tuple[str, int, int]:
    if _is_bed_row(fields):
        try:
            return fields[0], int(fields[1]), int(fields[2])
        except ValueError:
            pass
    if not fields:
        raise ValueError("empty selector row")
    match = _region_match(fields[0])
    if match is None:
        raise ValueError("{}:{} is not a region".format(selector, line_number))
    chrom, start, end = match.groups()
    return chrom, int(start), int(end)


def _text_score(fields: Sequence[str]) -> float:
    if len(fields) < 2:
        return 1.0
    for value in fields[1:]:
        try:
            return float(value)
        except ValueError:
            continue
    return 1.0


def _split_selector_fields(line: str) -> List[str]:
    stripped = line.strip()
    if "\t" in stripped:
        return stripped.split("\t")
    return stripped.split()


def _region_match(value: str):
    delimiters = re.escape(configured_peak_txt_delimiter())
    pattern = re.compile(
        r"^(.+?)[{}](\d+)[{}](\d+)$".format(delimiters, delimiters)
    )
    return pattern.match(value)


def _column_score(
    fields: Sequence[str],
    column: int,
    selector: PathLike,
    line_number: int,
) -> float:
    try:
        return float(fields[column - 1])
    except (IndexError, ValueError) as exc:
        raise ValueError(
            "{}:{} has no numeric score in configured column {}.".format(
                selector, line_number, column
            )
        ) from exc


def _sort_interval_map(intervals_by_chrom: Mapping[str, List[Tuple]]) -> Mapping:
    return {
        chrom: sorted(intervals, key=lambda item: (item[0], item[1]))
        for chrom, intervals in intervals_by_chrom.items()
    }


def _read_ranked_ids(selector: PathLike) -> Mapping[str, int]:
    ids = {}  # type: Dict[str, int]
    for identifier in _iter_selector_tokens(selector):
        ids.setdefault(identifier, len(ids))
    return ids


def _read_expression_table(selector: PathLike) -> Mapping[str, float]:
    expression = {}  # type: Dict[str, float]
    with Path(selector).expanduser().open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.strip().split()
            if len(fields) < 2:
                raise ValueError(
                    "{}:{} must contain isoform ID and expression.".format(
                        selector, line_number
                    )
                )
            try:
                value = float(fields[1])
            except ValueError as exc:
                if line_number == 1:
                    continue
                raise ValueError(
                    "{}:{} has non-numeric expression.".format(selector, line_number)
                ) from exc
            expression[fields[0]] = max(value, expression.get(fields[0], value))
    return expression


def _iter_selector_tokens(selector: PathLike) -> Iterable[str]:
    with Path(selector).expanduser().open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            for token in line.replace(",", " ").split():
                yield token


def _match_identifier(
    record: _IsoformRecord,
    ranked_ids: Mapping[str, int],
    inclusive: bool,
) -> Optional[int]:
    matches = [
        ranked_ids[identifier]
        for identifier in _record_identifiers(record, inclusive)
        if identifier in ranked_ids
    ]
    if matches:
        return min(matches)
    return None


def _record_identifiers(
    record: _IsoformRecord,
    inclusive: bool,
) -> Tuple[str, ...]:
    if not inclusive:
        return (record.transcript_id,)
    values = [record.transcript_id]
    if record.transcript_id_base and record.transcript_id_base not in values:
        values.append(record.transcript_id_base)
    stripped = _strip_version(record.transcript_id)
    if stripped and stripped not in values:
        values.append(stripped)
    return tuple(values)


def _merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not intervals:
        return []
    merged = []
    current_start, current_end = sorted(intervals)[0]
    for start, end in sorted(intervals)[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
            continue
        merged.append((current_start, current_end))
        current_start, current_end = start, end
    merged.append((current_start, current_end))
    return merged


def _normalize_method(method: str) -> str:
    aliases = {
        "longcol5": "longcol5",
        "long-col5": "longcol5",
        "long_col5": "longcol5",
        "long": "long",
        "length": "long",
        "peak": "peak",
        "peaks": "peak",
        "isoid": "isoID",
        "iso-id": "isoID",
        "iso_id": "isoID",
        "isoID": "isoID",
        "isoexp": "isoexp",
        "iso-exp": "isoexp",
        "iso_exp": "isoexp",
        "expression": "isoexp",
        "perover": "perover",
        "percent-overlap": "perover",
        "percent_overlap": "perover",
        "overlap": "perover",
    }
    try:
        return aliases[method]
    except KeyError:
        lowered = method.lower()
        if lowered in aliases:
            return aliases[lowered]
    raise ValueError(
        "Unsupported dedup/filter method {!r}. Supported: {}.".format(
            method, ", ".join(DEDUP_METHODS)
        )
    )


def _normalize_gene_key(gene_key: str) -> str:
    aliases = {
        "symbol": "symbol",
        "gene_symbol": "symbol",
        "gene-symbol": "symbol",
        "name": "symbol",
        "ensid": "ensid",
        "ens": "ensid",
        "ensembl": "ensid",
        "gene_id": "ensid",
        "gene-id": "ensid",
    }
    try:
        return aliases[gene_key.lower()]
    except KeyError as exc:
        raise ValueError("gene_key must be 'symbol' or 'ensid'.") from exc


def _parse_bp(value: Union[int, str]) -> int:
    if isinstance(value, int):
        return value
    text = value.strip().lower()
    multiplier = 1
    if text.endswith("mb"):
        multiplier = 1000000
        text = text[:-2]
    elif text.endswith("m"):
        multiplier = 1000000
        text = text[:-1]
    elif text.endswith("kb"):
        multiplier = 1000
        text = text[:-2]
    elif text.endswith("k"):
        multiplier = 1000
        text = text[:-1]
    elif text.endswith("bp"):
        text = text[:-2]
    parsed = int(float(text) * multiplier)
    if parsed <= 0:
        raise ValueError("promoter_bp must be positive.")
    return parsed


def _field_or_empty(fields: Sequence[str], index: int) -> str:
    if index < len(fields):
        return fields[index]
    return ""


def _strip_version(identifier: str) -> str:
    return identifier.split(".", 1)[0] if identifier else ""


# Short public names; retain the historical GENCODE names for compatibility.
dedup_bed = dedup_gencode_bed
filter_bed = filter_gencode_bed
