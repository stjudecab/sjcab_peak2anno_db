"""Generate merged GENCODE region BED classes."""

from __future__ import annotations

import json
import os
import re
import shutil
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Tuple, Union

from ._download import ProgressCallback, report_progress
from ._derive import _parse_bed_fields, _site_interval
from ._gencode import (
    _ENSEMBL_SPECIES,
    _open_text,
    _parse_attributes,
    _download_ucsc_primary_chromosomes,
    _ucsc_gtf_builds,
    convert_gencode_gtf_to_bed,
    download_gencode_gtf,
    gencode_bed_filename,
)
from ._registry import data_root, user_data_dir

PathLike = Union[str, os.PathLike]

TSS_FLANK_REGION_TYPES = ("dis5", "dis3", "promoter.up", "promoter.down")
GENCODE_REGION_TYPES = ("exon", "intron", "tes", "intergenic")
REGION_TYPES = TSS_FLANK_REGION_TYPES + GENCODE_REGION_TYPES
GENCODE_FEATURE_LIST_ORDER = (
    "promoter.up",
    "promoter.down",
    "exon",
    "intron",
    "tes",
    "dis5",
    "dis3",
    "intergenic",
)


def write_tss_flank_region_unions(
    gene_bed: PathLike,
    output_dir: PathLike,
    promoter_bp: Union[int, str] = 2000,
    distal_bp: Union[int, str] = 50000,
    prefix: Optional[str] = None,
    split_tss: bool = True,
) -> Mapping[str, Path]:
    """Write merged old-style promoter/distal BED classes from a gene BED file.

    The generated files match the legacy CAB/``annotate_prep.sh`` classes:
    promoter and ``dis5`` are strand-aware TSS flanks, while ``dis3`` is a
    strand-aware TES downstream flank. By default, ``split_tss=True`` models the
    1 bp TSS/TES interval used by BEDTools ``flankBed``, so the site base is
    excluded before merging.
    """

    promoter = _parse_bp(promoter_bp)
    distal = _parse_bp(distal_bp)
    if promoter <= 0:
        raise ValueError("promoter_bp must be positive.")
    if distal <= promoter:
        raise ValueError("distal_bp must be greater than promoter_bp.")

    label = prefix or _format_bp(promoter)
    intervals, chrom_order = _collect_tss_flank_intervals(
        gene_bed,
        promoter,
        distal,
        split_tss=split_tss,
    )

    target_dir = Path(output_dir).expanduser()
    target_dir.mkdir(parents=True, exist_ok=True)
    outputs = {}  # type: Dict[str, Path]
    for region_type in TSS_FLANK_REGION_TYPES:
        output_path = target_dir / "{}.{}.bed".format(label, region_type)
        _write_merged_intervals(
            intervals[region_type],
            chrom_order,
            output_path,
            merge_gap=2,
        )
        outputs[region_type] = output_path
    return outputs


def write_gencode_region_unions(
    gtf_path: PathLike,
    output_dir: PathLike,
    tes_bp: Union[int, str] = 2000,
    prefix: Optional[str] = None,
) -> Mapping[str, Path]:
    """Write merged exon, isoform intron, TES-window, and intergenic BED classes.

    ``intron.bed`` is built from the gaps between adjacent exon records for
    each transcript isoform. ``tes.bed`` covers TES - ``tes_bp`` through TES +
    ``tes_bp`` and includes the TES base. ``intergenic.bed`` is the complement
    of merged gene intervals using GTF ``##sequence-region`` chromosome lengths
    when present; otherwise it uses the observed annotation span for each
    chromosome.
    """

    flank = _parse_bp(tes_bp)
    if flank <= 0:
        raise ValueError("tes_bp must be positive.")

    parsed = _collect_gencode_regions(gtf_path, flank)
    target_dir = Path(output_dir).expanduser()
    target_dir.mkdir(parents=True, exist_ok=True)

    outputs = {}  # type: Dict[str, Path]
    for region_type in GENCODE_REGION_TYPES:
        if prefix:
            filename = "{}.{}.bed".format(prefix, region_type)
        else:
            filename = "{}.bed".format(region_type)
        output_path = target_dir / filename
        _write_merged_intervals(
            parsed.intervals[region_type],
            parsed.chrom_order,
            output_path,
        )
        outputs[region_type] = output_path
    return outputs


def download_gencode_feature(
    species: str,
    version: str,
    output_dir: PathLike = ".",
    gtf_path: Optional[PathLike] = None,
    gtf_url: Optional[str] = None,
    gene_bed: Optional[PathLike] = None,
    promoter_bp: Union[int, str] = 2000,
    distal_bp: Union[int, str] = 50000,
    prefix: Optional[str] = None,
    split_tss: bool = True,
    include_type: bool = True,
    tes_bp: Union[int, str] = 2000,
    overwrite: bool = True,
    log_data_dir: Optional[PathLike] = None,
    progress: Optional[ProgressCallback] = None,
    source: str = "auto",
    cache_dir: Optional[PathLike] = None,
    ucsc_annotation: str = "ens",
    clean_cache: bool = False,
    *,
    gene_bed_output: Optional[PathLike] = None,
) -> Mapping[str, Path]:
    """Download/convert GENCODE GTF and write old-style feature BED classes."""

    target_dir = Path(output_dir).expanduser()
    target_dir.mkdir(parents=True, exist_ok=True)
    report_progress(
        progress,
        "download-gencode-feature {} {}: started".format(species, version),
    )

    report_progress(progress, "download-gencode-feature: preparing GTF")
    gtf_for_regions = _resolve_gtf_for_regions(
        species,
        version,
        target_dir,
        gtf_path=gtf_path,
        gtf_url=gtf_url,
        overwrite=overwrite,
        log_data_dir=log_data_dir,
        progress=progress,
        source=source,
        cache_dir=cache_dir,
        ucsc_annotation=ucsc_annotation,
    )
    report_progress(progress, "download-gencode-feature: GTF ready")

    if gene_bed is None:
        if gene_bed_output is None:
            gene_bed = target_dir / gencode_bed_filename(
                species, version, include_type=include_type
            )
        else:
            gene_bed = Path(gene_bed_output).expanduser()
        if overwrite or not Path(gene_bed).exists():
            report_progress(progress, "download-gencode-feature: converting gene BED")
            convert_gencode_gtf_to_bed(
                gtf_for_regions,
                gene_bed,
                include_type=include_type,
            )
            report_progress(progress, "download-gencode-feature: gene BED done")
    else:
        gene_bed = Path(gene_bed).expanduser()
        report_progress(
            progress,
            "download-gencode-feature: using existing gene BED {}".format(gene_bed),
        )

    label = gencode_feature_prefix(promoter_bp=promoter_bp, prefix=prefix)
    report_progress(
        progress,
        "download-gencode-feature: writing legacy feature BED files",
    )
    outputs = write_legacy_gencode_feature_unions(
        gene_bed,
        gtf_for_regions,
        target_dir,
        promoter_bp=promoter_bp,
        distal_bp=distal_bp,
        tes_bp=tes_bp,
        prefix=label,
        split_tss=split_tss,
        chrom_sizes=_resolve_species_chrom_sizes(species),
    )
    report_progress(progress, "download-gencode-feature: legacy BED files done")
    outputs["list"] = write_gencode_feature_list(outputs, target_dir, label)
    if clean_cache and gtf_path is None:
        cached_path = Path(gtf_for_regions)
        cache_root = user_data_dir(cache_dir) / "cachegtf"
        if cached_path.parent == cache_root and cached_path.exists():
            cached_path.unlink()
    report_progress(progress, "download-gencode-feature: feature list done")
    outputs["gene_bed"] = Path(gene_bed).expanduser()
    report_progress(
        progress,
        "download-gencode-feature {} {}: done".format(species, version),
    )
    return outputs


download_feature = download_gencode_feature


def download_gencode_tss_flank_region_unions(
    species: str,
    version: str,
    output_dir: PathLike = ".",
    gtf_path: Optional[PathLike] = None,
    gtf_url: Optional[str] = None,
    gene_bed: Optional[PathLike] = None,
    promoter_bp: Union[int, str] = 2000,
    distal_bp: Union[int, str] = 50000,
    prefix: Optional[str] = None,
    split_tss: bool = True,
    include_type: bool = True,
    tes_bp: Union[int, str] = 2000,
    overwrite: bool = True,
    log_data_dir: Optional[PathLike] = None,
    progress: Optional[ProgressCallback] = None,
) -> Mapping[str, Path]:
    """Backward-compatible alias for :func:`download_gencode_feature`."""

    return download_gencode_feature(
        species,
        version,
        output_dir,
        gtf_path=gtf_path,
        gtf_url=gtf_url,
        gene_bed=gene_bed,
        promoter_bp=promoter_bp,
        distal_bp=distal_bp,
        prefix=prefix,
        split_tss=split_tss,
        include_type=include_type,
        tes_bp=tes_bp,
        overwrite=overwrite,
        log_data_dir=log_data_dir,
        progress=progress,
    )


def write_gencode_feature_list(
    outputs: Mapping[str, Path],
    output_dir: PathLike,
    prefix: str,
) -> Path:
    """Write an ordered list of generated GENCODE feature BED files."""

    target_dir = Path(output_dir).expanduser()
    list_path = target_dir / "order.lst"
    tmp_path = list_path.with_name(list_path.name + ".tmp")
    try:
        with tmp_path.open("w", encoding="utf-8") as handle:
            for region_type in GENCODE_FEATURE_LIST_ORDER:
                handle.write("{}\n".format(outputs[region_type].name))
        tmp_path.replace(list_path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()
    return list_path


def gencode_feature_prefix(
    promoter_bp: Union[int, str] = 2000,
    prefix: Optional[str] = None,
) -> str:
    """Return the output prefix for one GENCODE feature set."""

    return prefix or _format_bp(_parse_bp(promoter_bp))


def write_legacy_gencode_feature_unions(
    gene_bed: PathLike,
    gtf_path: PathLike,
    output_dir: PathLike,
    promoter_bp: Union[int, str] = 2000,
    distal_bp: Union[int, str] = 50000,
    tes_bp: Union[int, str] = 2000,
    prefix: Optional[str] = None,
    split_tss: bool = True,
    chrom_sizes: Optional[PathLike] = None,
) -> Mapping[str, Path]:
    """Write legacy CAB/``annotate_prep.sh`` feature BED classes.

    This mirrors the old preparation script without requiring BEDTools:
    promoter and ``dis5`` are strand-aware TSS flanks, ``dis3`` is a
    strand-aware downstream TES flank, ``tes`` is the merged two-sided TES
    flank, introns overlapping more than ten merged promoter-up intervals are
    removed, and intergenic is the complement of all generated feature classes.
    """

    promoter = _parse_bp(promoter_bp)
    distal = _parse_bp(distal_bp)
    tes = _parse_bp(tes_bp)
    if promoter <= 0:
        raise ValueError("promoter_bp must be positive.")
    if distal <= promoter:
        raise ValueError("distal_bp must be greater than promoter_bp.")
    if tes <= 0:
        raise ValueError("tes_bp must be positive.")

    label = prefix or _format_bp(promoter)
    parsed = _collect_legacy_gencode_feature_regions(
        gene_bed,
        gtf_path,
        promoter=promoter,
        distal=distal,
        tes=tes,
        split_tss=split_tss,
        fallback_chrom_lengths=_read_chrom_sizes(chrom_sizes),
    )

    target_dir = Path(output_dir).expanduser()
    target_dir.mkdir(parents=True, exist_ok=True)
    outputs = {}  # type: Dict[str, Path]
    for region_type in GENCODE_FEATURE_LIST_ORDER + ("promoter",):
        output_path = target_dir / "{}.{}.bed".format(label, region_type)
        _write_intervals(
            parsed.intervals[region_type],
            parsed.chrom_order,
            output_path,
        )
        outputs[region_type] = output_path
    return outputs


class _ParsedGencodeRegions:
    def __init__(
        self,
        intervals: Dict[str, Dict[str, List[Tuple[int, int]]]],
        chrom_order: Tuple[str, ...],
    ) -> None:
        self.intervals = intervals
        self.chrom_order = chrom_order


def _resolve_gtf_for_regions(
    species: str,
    version: str,
    target_dir: Path,
    gtf_path: Optional[PathLike],
    gtf_url: Optional[str],
    overwrite: bool,
    log_data_dir: Optional[PathLike],
    progress: Optional[ProgressCallback],
    source: str,
    cache_dir: Optional[PathLike],
    ucsc_annotation: str,
) -> Path:
    if gtf_path is not None:
        return Path(gtf_path).expanduser()
    return download_gencode_gtf(
        species,
        version,
        target_dir,
        gtf_url=gtf_url,
        overwrite=overwrite,
        log_data_dir=log_data_dir,
        progress=progress,
        source=source,
        cache_dir=cache_dir,
        ucsc_annotation=ucsc_annotation,
    )


def _resolve_species_chrom_sizes(species: str) -> Optional[Path]:
    cache_dir = user_data_dir() / "sizes"
    cache_dir.mkdir(parents=True, exist_ok=True)
    raw_path = cache_dir / "{}.sizes".format(species)
    clean_path = cache_dir / "{}.sizes.clean".format(species)

    packaged_clean = _cache_packaged_sizes(species, raw_path)
    if not raw_path.exists():
        try:
            if species in _ucsc_gtf_builds():
                _download_ucsc_sizes(species, raw_path)
            else:
                _download_ensembl_sizes(species, raw_path)
        except (OSError, ValueError, KeyError, json.JSONDecodeError):
            pass
    candidates = (
        Path("~/data/{}.sizes".format(species)).expanduser(),
        Path("~/data/{}/{}.sizes".format(species, species)).expanduser(),
    )
    if not raw_path.exists():
        for candidate in candidates:
            if candidate.exists():
                shutil.copyfile(candidate, raw_path)
                break
    if not raw_path.exists():
        return None

    if not clean_path.exists():
        if packaged_clean is not None and packaged_clean.exists():
            shutil.copyfile(packaged_clean, clean_path)
        elif species in _ucsc_gtf_builds():
            _write_clean_chrom_sizes(
                raw_path, clean_path, _download_ucsc_primary_chromosomes(species)
            )
        else:
            _write_clean_chrom_sizes(
                raw_path, clean_path, _download_ensembl_primary_chromosomes(species)
            )
    return clean_path


def _cache_packaged_sizes(species: str, destination: Path) -> Optional[Path]:
    species_key = species.lower().replace(" ", "_")
    names = [species]
    metadata = _ENSEMBL_SPECIES.get(species_key)
    if metadata is not None:
        names.insert(0, metadata[0])
    for name in names:
        source = data_root() / "sizes" / "{}.sizes".format(name)
        if source.exists():
            shutil.copyfile(source, destination)
            clean_source = source.with_name(source.name + ".clean")
            return clean_source if clean_source.exists() else None
    return None


def _download_ucsc_sizes(species: str, destination: Path) -> None:
    url = "https://hgdownload.soe.ucsc.edu/goldenPath/{0}/bigZips/{0}.chrom.sizes".format(
        urllib.parse.quote(species, safe="")
    )
    _write_sizes_from_text(url, destination)


def _download_ensembl_sizes(species: str, destination: Path) -> None:
    url = "https://rest.ensembl.org/info/assembly/{}?content-type=application/json".format(
        urllib.parse.quote(species, safe="")
    )
    payload = _fetch_json(url)
    rows = payload.get("top_level_region", [])
    _write_sizes(
        destination,
        ((row["name"], row["length"]) for row in rows if "name" in row),
    )


def _download_ensembl_primary_chromosomes(species: str):
    url = "https://rest.ensembl.org/info/assembly/{}?content-type=application/json".format(
        urllib.parse.quote(species, safe="")
    )
    payload = _fetch_json(url)
    return {
        row["name"]
        for row in payload.get("top_level_region", [])
        if row.get("coord_system") == "chromosome" and "name" in row
    }


def _fetch_json(url: str):
    request = urllib.request.Request(
        url, headers={"User-Agent": "sjcab_peak2anno_db"}
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.loads(response.read().decode("utf-8"))


def _write_sizes_from_text(url: str, destination: Path) -> None:
    request = urllib.request.Request(
        url, headers={"User-Agent": "sjcab_peak2anno_db"}
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        rows = response.read().decode("utf-8").splitlines()
    _write_sizes(
        destination,
        (line.split()[:2] for line in rows if len(line.split()) >= 2),
    )


def _write_sizes(destination: Path, rows: Iterable[Tuple[str, int]]) -> None:
    temporary = destination.with_name(destination.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        for name, size in rows:
            handle.write("{}\t{}\n".format(name, int(size)))
    temporary.replace(destination)


def _write_clean_chrom_sizes(
    source: Path, destination: Path, primary_chromosomes=None
) -> None:
    temporary = destination.with_name(destination.name + ".tmp")
    with source.open("r", encoding="utf-8") as source_handle, temporary.open(
        "w", encoding="utf-8"
    ) as destination_handle:
        for line in source_handle:
            fields = line.rstrip("\n").split()
            if len(fields) >= 2 and (
                (not primary_chromosomes and _is_primary_chromosome(fields[0]))
                or (primary_chromosomes and fields[0] in primary_chromosomes)
            ):
                destination_handle.write("{}\t{}\n".format(fields[0], fields[1]))
    temporary.replace(destination)


def _is_primary_chromosome(name: str) -> bool:
    match = re.fullmatch(r"(?:chr)?(\d+|X|Y|M|MT)", name)
    if match is None:
        return False
    suffix = match.group(1)
    return suffix in {"X", "Y", "M", "MT"} or 1 <= int(suffix) <= 22


def _read_chrom_sizes(
    chrom_sizes: Optional[PathLike],
) -> Mapping[str, int]:
    if chrom_sizes is None:
        return {}
    sizes_path = Path(chrom_sizes).expanduser()
    lengths = {}  # type: Dict[str, int]
    with sizes_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split()
            if len(fields) < 2:
                raise ValueError(
                    "{}:{} is not a valid chromosome sizes row".format(
                        sizes_path, line_number
                    )
                )
            try:
                lengths[fields[0]] = int(fields[1])
            except ValueError as exc:
                raise ValueError(
                    "{}:{} has a non-integer chromosome size".format(
                        sizes_path, line_number
                    )
                ) from exc
    return lengths


def _collect_gencode_regions(
    gtf_path: PathLike,
    tes_bp: int,
) -> _ParsedGencodeRegions:
    intervals = {
        region_type: {} for region_type in GENCODE_REGION_TYPES
    }  # type: Dict[str, Dict[str, List[Tuple[int, int]]]]
    chrom_seen = {}  # type: Dict[str, None]
    chrom_lengths = {}  # type: Dict[str, int]
    chrom_max_end = {}  # type: Dict[str, int]
    gene_intervals = {}  # type: Dict[str, List[Tuple[int, int]]]
    transcripts = {}  # type: Dict[str, Tuple[str, int, int, str]]
    transcript_exons = {}  # type: Dict[str, List[Tuple[str, int, int]]]

    with _open_text(gtf_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            if line.startswith("##sequence-region"):
                _record_sequence_region(line, chrom_seen, chrom_lengths)
                continue
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                raise ValueError(
                    "{}:{} is not a valid GTF row".format(gtf_path, line_number)
                )

            feature = fields[2]
            if feature not in {"gene", "transcript", "exon"}:
                continue

            chrom = fields[0]
            start = int(fields[3]) - 1
            end = int(fields[4])
            strand = fields[6]
            attrs = _parse_attributes(fields[8])

            if chrom not in chrom_seen:
                chrom_seen[chrom] = None
            chrom_max_end[chrom] = max(chrom_max_end.get(chrom, 0), end)

            if feature == "gene":
                gene_intervals.setdefault(chrom, []).append((start, end))
                continue

            transcript_id = attrs.get("transcript_id")
            if not transcript_id:
                continue

            if feature == "transcript":
                transcripts[transcript_id] = (chrom, start, end, strand)
                continue

            intervals["exon"].setdefault(chrom, []).append((start, end))
            transcript_exons.setdefault(transcript_id, []).append((chrom, start, end))
            existing = transcripts.get(transcript_id)
            if existing is None:
                transcripts[transcript_id] = (chrom, start, end, strand)
            else:
                _chrom, previous_start, previous_end, previous_strand = existing
                transcripts[transcript_id] = (
                    chrom,
                    min(previous_start, start),
                    max(previous_end, end),
                    previous_strand,
                )

    if not any(gene_intervals.values()):
        for chrom, start, end, _strand in transcripts.values():
            gene_intervals.setdefault(chrom, []).append((start, end))

    _collect_tes_windows(transcripts, intervals["tes"], tes_bp)
    _collect_isoform_introns(transcript_exons, intervals["intron"])
    _collect_intergenic(
        gene_intervals,
        chrom_lengths,
        chrom_max_end,
        tuple(chrom_seen.keys()),
        intervals["intergenic"],
    )

    return _ParsedGencodeRegions(intervals, tuple(chrom_seen.keys()))


def _collect_legacy_gencode_feature_regions(
    gene_bed: PathLike,
    gtf_path: PathLike,
    promoter: int,
    distal: int,
    tes: int,
    split_tss: bool,
    fallback_chrom_lengths: Mapping[str, int],
) -> _ParsedGencodeRegions:
    intervals = {
        region_type: {}
        for region_type in GENCODE_FEATURE_LIST_ORDER + ("promoter",)
    }  # type: Dict[str, Dict[str, List[Tuple[int, int]]]]
    intervals["_dis3_inner"] = {}
    chrom_seen = {}  # type: Dict[str, None]
    chrom_lengths = {}  # type: Dict[str, int]
    chrom_max_end = {}  # type: Dict[str, int]
    transcript_exons = {}  # type: Dict[str, List[Tuple[str, int, int]]]

    _collect_legacy_gtf_regions(
        gtf_path,
        intervals,
        transcript_exons,
        chrom_seen,
        chrom_lengths,
        chrom_max_end,
    )
    if not chrom_lengths:
        for chrom, length in fallback_chrom_lengths.items():
            if chrom not in chrom_seen:
                chrom_seen[chrom] = None
            chrom_lengths[chrom] = length
    _collect_isoform_introns(transcript_exons, intervals["intron"])
    _collect_legacy_gene_bed_regions(
        gene_bed,
        intervals,
        chrom_seen,
        chrom_lengths,
        chrom_max_end,
        promoter=promoter,
        distal=distal,
        tes=tes,
        split_tss=split_tss,
    )

    chrom_order = tuple(chrom_seen.keys())
    merged = {}  # type: Dict[str, Dict[str, List[Tuple[int, int]]]]
    dis5_intervals = _subtract_interval_map(
        intervals["dis5"],
        intervals["promoter.up"],
        chrom_order,
    )
    dis3_intervals = _subtract_interval_map(
        intervals["dis3"],
        intervals["_dis3_inner"],
        chrom_order,
    )
    for region_type in ("promoter.up", "promoter.down", "exon", "tes"):
        merged[region_type] = _merge_interval_map(
            intervals[region_type],
            chrom_order,
            max_gap=2,
        )
    merged["dis5"] = _merge_interval_map(dis5_intervals, chrom_order, max_gap=2)
    merged["dis3"] = _merge_interval_map(dis3_intervals, chrom_order, max_gap=2)

    promoter_source = _combine_interval_maps(
        intervals["promoter.up"],
        intervals["promoter.down"],
    )
    merged["promoter"] = _merge_interval_map(
        promoter_source,
        chrom_order,
        max_gap=2,
    )

    merged_introns = _merge_interval_map(
        intervals["intron"],
        chrom_order,
        max_gap=2,
    )
    merged["intron"] = _filter_introns_by_promoter_up(
        merged_introns,
        merged["promoter.up"],
        max_overlaps=10,
    )

    feature_union = _combine_interval_maps(
        merged["promoter"],
        merged["intron"],
        merged["exon"],
        merged["dis3"],
        merged["dis5"],
        merged["tes"],
    )
    merged_union = _merge_interval_map(feature_union, chrom_order, max_gap=2)
    merged["intergenic"] = _complement_merged_intervals(
        merged_union,
        chrom_lengths,
        chrom_max_end,
        chrom_order,
    )

    return _ParsedGencodeRegions(merged, chrom_order)


def _collect_legacy_gtf_regions(
    gtf_path: PathLike,
    intervals: Dict[str, Dict[str, List[Tuple[int, int]]]],
    transcript_exons: Dict[str, List[Tuple[str, int, int]]],
    chrom_seen: Dict[str, None],
    chrom_lengths: Dict[str, int],
    chrom_max_end: Dict[str, int],
) -> None:
    with _open_text(gtf_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            if line.startswith("##sequence-region"):
                _record_sequence_region(line, chrom_seen, chrom_lengths)
                continue
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                raise ValueError(
                    "{}:{} is not a valid GTF row".format(gtf_path, line_number)
                )

            feature = fields[2]
            if feature != "exon":
                continue

            chrom = fields[0]
            start = int(fields[3]) - 1
            end = int(fields[4])
            attrs = _parse_attributes(fields[8])
            transcript_id = attrs.get("transcript_id")
            if not transcript_id:
                continue

            if chrom not in chrom_seen:
                chrom_seen[chrom] = None
            _record_interval_max_end(chrom_max_end, chrom, end)
            intervals["exon"].setdefault(chrom, []).append((start, end))
            transcript_exons.setdefault(transcript_id, []).append((chrom, start, end))


def _collect_legacy_gene_bed_regions(
    gene_bed: PathLike,
    intervals: Dict[str, Dict[str, List[Tuple[int, int]]]],
    chrom_seen: Dict[str, None],
    chrom_lengths: Mapping[str, int],
    chrom_max_end: Dict[str, int],
    promoter: int,
    distal: int,
    tes: int,
    split_tss: bool,
) -> None:
    with Path(gene_bed).expanduser().open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = _parse_bed_fields(line, line_number)
            chrom = fields[0]
            strand = fields[5]
            chrom_end = chrom_lengths.get(chrom)
            if chrom not in chrom_seen:
                chrom_seen[chrom] = None

            tss_start, tss_end = _site_span(fields, "tss", split_tss)
            tes_start, tes_end = _site_span(fields, "tes", split_tss)

            promoter_up = _upstream_site_flank(
                tss_start,
                tss_end,
                strand,
                promoter,
                chrom_end,
            )
            promoter_down = _downstream_site_flank(
                tss_start,
                tss_end,
                strand,
                promoter,
                chrom_end,
            )
            _add_interval(intervals["promoter.up"], chrom, promoter_up)
            _add_interval(intervals["promoter.down"], chrom, promoter_down)

            distal_up = _upstream_site_flank(
                tss_start,
                tss_end,
                strand,
                distal,
                chrom_end,
            )
            _add_interval(intervals["dis5"], chrom, distal_up)

            distal_down = _downstream_site_flank(
                tes_start,
                tes_end,
                strand,
                distal,
                chrom_end,
            )
            promoter_down_from_tes = _downstream_site_flank(
                tes_start,
                tes_end,
                strand,
                promoter,
                chrom_end,
            )
            _add_interval(intervals["dis3"], chrom, distal_down)
            _add_interval(intervals["_dis3_inner"], chrom, promoter_down_from_tes)

            _add_interval(
                intervals["tes"],
                chrom,
                _upstream_site_flank(tes_start, tes_end, strand, tes, chrom_end),
            )
            _add_interval(
                intervals["tes"],
                chrom,
                _downstream_site_flank(tes_start, tes_end, strand, tes, chrom_end),
            )

            for region_type in (
                "promoter.up",
                "promoter.down",
                "dis5",
                "dis3",
                "tes",
            ):
                for _start, end in intervals[region_type].get(chrom, [])[-2:]:
                    _record_interval_max_end(chrom_max_end, chrom, end)


def _record_sequence_region(
    line: str,
    chrom_seen: Dict[str, None],
    chrom_lengths: Dict[str, int],
) -> None:
    fields = line.strip().split()
    if len(fields) < 4:
        return
    chrom = fields[1]
    try:
        end = int(fields[3])
    except ValueError:
        return
    if chrom not in chrom_seen:
        chrom_seen[chrom] = None
    chrom_lengths[chrom] = max(chrom_lengths.get(chrom, 0), end)


def _collect_tes_windows(
    transcripts: Mapping[str, Tuple[str, int, int, str]],
    tes_intervals: Dict[str, List[Tuple[int, int]]],
    flank: int,
) -> None:
    for chrom, start, end, strand in transcripts.values():
        point = start if strand == "-" else end - 1
        point = max(point, 0)
        tes_intervals.setdefault(chrom, []).append(
            (max(0, point - flank), point + flank + 1)
        )


def _collect_isoform_introns(
    transcript_exons: Mapping[str, List[Tuple[str, int, int]]],
    intron_intervals: Dict[str, List[Tuple[int, int]]],
) -> None:
    for exons in transcript_exons.values():
        exons_by_chrom = {}  # type: Dict[str, List[Tuple[int, int]]]
        for chrom, start, end in exons:
            exons_by_chrom.setdefault(chrom, []).append((start, end))

        for chrom, chrom_exons in exons_by_chrom.items():
            merged_exons = _merge_intervals(chrom_exons)
            for previous, current in zip(merged_exons, merged_exons[1:]):
                intron_start = previous[1]
                intron_end = current[0]
                if intron_end > intron_start:
                    intron_intervals.setdefault(chrom, []).append(
                        (intron_start, intron_end)
                    )


def _collect_intergenic(
    gene_intervals: Mapping[str, List[Tuple[int, int]]],
    chrom_lengths: Mapping[str, int],
    chrom_max_end: Mapping[str, int],
    chrom_order: Tuple[str, ...],
    intergenic_intervals: Dict[str, List[Tuple[int, int]]],
) -> None:
    for chrom in chrom_order:
        chrom_end = chrom_lengths.get(chrom, chrom_max_end.get(chrom, 0))
        if chrom_end <= 0:
            continue
        cursor = 0
        for start, end in _merge_intervals(gene_intervals.get(chrom, [])):
            if start > cursor:
                intergenic_intervals.setdefault(chrom, []).append((cursor, start))
            cursor = max(cursor, end)
        if cursor < chrom_end:
            intergenic_intervals.setdefault(chrom, []).append((cursor, chrom_end))


def _complement_merged_intervals(
    occupied_by_chrom: Mapping[str, List[Tuple[int, int]]],
    chrom_lengths: Mapping[str, int],
    chrom_max_end: Mapping[str, int],
    chrom_order: Tuple[str, ...],
) -> Dict[str, List[Tuple[int, int]]]:
    complements = {}  # type: Dict[str, List[Tuple[int, int]]]
    for chrom in chrom_order:
        chrom_end = chrom_lengths.get(chrom, chrom_max_end.get(chrom, 0))
        if chrom_end <= 0:
            continue
        cursor = 0
        for start, end in occupied_by_chrom.get(chrom, []):
            if start > cursor:
                complements.setdefault(chrom, []).append((cursor, start))
            cursor = max(cursor, end)
        if cursor < chrom_end:
            complements.setdefault(chrom, []).append((cursor, chrom_end))
    return complements


def _collect_tss_flank_intervals(
    gene_bed: PathLike,
    promoter: int,
    distal: int,
    split_tss: bool,
) -> Tuple[Dict[str, Dict[str, List[Tuple[int, int]]]], Tuple[str, ...]]:
    intervals = {
        region_type: {} for region_type in TSS_FLANK_REGION_TYPES
    }  # type: Dict[str, Dict[str, List[Tuple[int, int]]]]
    dis3_inner = {}  # type: Dict[str, List[Tuple[int, int]]]
    chrom_seen = {}  # type: Dict[str, None]

    with Path(gene_bed).expanduser().open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = _parse_bed_fields(line, line_number)
            chrom = fields[0]
            strand = fields[5]
            if chrom not in chrom_seen:
                chrom_seen[chrom] = None

            tss_start, tss_end = _site_span(fields, "tss", split_tss)
            tes_start, tes_end = _site_span(fields, "tes", split_tss)
            promoter_up = _upstream_site_flank(
                tss_start,
                tss_end,
                strand,
                promoter,
                chrom_end=None,
            )
            promoter_down = _downstream_site_flank(
                tss_start,
                tss_end,
                strand,
                promoter,
                chrom_end=None,
            )
            _add_interval(intervals["promoter.up"], chrom, promoter_up)
            _add_interval(intervals["promoter.down"], chrom, promoter_down)

            distal_up = _upstream_site_flank(
                tss_start,
                tss_end,
                strand,
                distal,
                chrom_end=None,
            )
            _add_interval(intervals["dis5"], chrom, distal_up)

            distal_down = _downstream_site_flank(
                tes_start,
                tes_end,
                strand,
                distal,
                chrom_end=None,
            )
            promoter_down_from_tes = _downstream_site_flank(
                tes_start,
                tes_end,
                strand,
                promoter,
                chrom_end=None,
            )
            _add_interval(intervals["dis3"], chrom, distal_down)
            _add_interval(dis3_inner, chrom, promoter_down_from_tes)

    chrom_order = tuple(chrom_seen.keys())
    intervals["dis5"] = _subtract_interval_map(
        intervals["dis5"],
        intervals["promoter.up"],
        chrom_order,
    )
    intervals["dis3"] = _subtract_interval_map(
        intervals["dis3"],
        dis3_inner,
        chrom_order,
    )
    return intervals, chrom_order


def _tss_flank_intervals(
    tss: int,
    promoter: int,
    distal: int,
    split_tss: bool,
) -> Iterable[Tuple[str, int, int]]:
    left_distal_start = max(0, tss - distal)
    left_promoter_start = max(0, tss - promoter)
    right_start = tss + 1 if split_tss else tss
    right_promoter_end = right_start + promoter
    right_distal_start = right_promoter_end
    right_distal_end = right_start + distal

    yield "dis5", left_distal_start, left_promoter_start
    yield "promoter.up", left_promoter_start, tss
    yield "promoter.down", right_start, right_promoter_end
    yield "dis3", right_distal_start, right_distal_end


def _site_span(
    fields: List[str],
    site: str,
    split_site: bool,
) -> Tuple[int, int]:
    start_text, end_text = _site_interval(fields, site)
    start = int(start_text)
    if split_site:
        return start, int(end_text)
    return start, start


def _upstream_site_flank(
    site_start: int,
    site_end: int,
    strand: str,
    bp: int,
    chrom_end: Optional[int],
) -> Tuple[int, int]:
    if strand == "-":
        return _clip_interval(site_end, site_end + bp, chrom_end)
    return _clip_interval(site_start - bp, site_start, chrom_end)


def _downstream_site_flank(
    site_start: int,
    site_end: int,
    strand: str,
    bp: int,
    chrom_end: Optional[int],
) -> Tuple[int, int]:
    if strand == "-":
        return _clip_interval(site_start - bp, site_start, chrom_end)
    return _clip_interval(site_end, site_end + bp, chrom_end)


def _clip_interval(
    start: int,
    end: int,
    chrom_end: Optional[int],
) -> Tuple[int, int]:
    start = max(0, start)
    end = max(0, end)
    if chrom_end is not None:
        start = min(start, chrom_end)
        end = min(end, chrom_end)
    return start, end


def _subtract_interval(
    interval: Tuple[int, int],
    subtract: Tuple[int, int],
) -> List[Tuple[int, int]]:
    start, end = interval
    subtract_start, subtract_end = subtract
    if end <= start:
        return []
    if subtract_end <= start or subtract_start >= end:
        return [(start, end)]
    pieces = []
    left = (start, max(start, min(end, subtract_start)))
    right = (min(end, max(start, subtract_end)), end)
    if left[1] > left[0]:
        pieces.append(left)
    if right[1] > right[0]:
        pieces.append(right)
    return pieces


def _subtract_interval_map(
    intervals_by_chrom: Mapping[str, List[Tuple[int, int]]],
    subtract_by_chrom: Mapping[str, List[Tuple[int, int]]],
    chrom_order: Tuple[str, ...],
) -> Dict[str, List[Tuple[int, int]]]:
    subtracted = {}  # type: Dict[str, List[Tuple[int, int]]]
    for chrom in chrom_order:
        intervals = sorted(intervals_by_chrom.get(chrom, []))
        subtractors = sorted(subtract_by_chrom.get(chrom, []))
        first_subtractor = 0
        chrom_pieces = []
        for start, end in intervals:
            if end <= start:
                continue
            while (
                first_subtractor < len(subtractors)
                and subtractors[first_subtractor][1] <= start
            ):
                first_subtractor += 1
            pieces = [(start, end)]
            index = first_subtractor
            while index < len(subtractors) and subtractors[index][0] < end:
                next_pieces = []
                for piece in pieces:
                    next_pieces.extend(_subtract_interval(piece, subtractors[index]))
                pieces = next_pieces
                if not pieces:
                    break
                index += 1
            chrom_pieces.extend(pieces)
        if chrom_pieces:
            subtracted[chrom] = chrom_pieces
    return subtracted


def _add_interval(
    intervals_by_chrom: Dict[str, List[Tuple[int, int]]],
    chrom: str,
    interval: Tuple[int, int],
) -> None:
    start, end = interval
    if end > start:
        intervals_by_chrom.setdefault(chrom, []).append(interval)


def _record_interval_max_end(
    chrom_max_end: Dict[str, int],
    chrom: str,
    end: int,
) -> None:
    chrom_max_end[chrom] = max(chrom_max_end.get(chrom, 0), end)


def _combine_interval_maps(
    *interval_maps: Mapping[str, List[Tuple[int, int]]]
) -> Dict[str, List[Tuple[int, int]]]:
    combined = {}  # type: Dict[str, List[Tuple[int, int]]]
    for interval_map in interval_maps:
        for chrom, intervals in interval_map.items():
            combined.setdefault(chrom, []).extend(intervals)
    return combined


def _merge_interval_map(
    intervals_by_chrom: Mapping[str, List[Tuple[int, int]]],
    chrom_order: Tuple[str, ...],
    max_gap: int = 0,
) -> Dict[str, List[Tuple[int, int]]]:
    merged = {}  # type: Dict[str, List[Tuple[int, int]]]
    for chrom in chrom_order:
        chrom_intervals = _merge_intervals(
            list(intervals_by_chrom.get(chrom, [])),
            max_gap=max_gap,
        )
        if chrom_intervals:
            merged[chrom] = chrom_intervals
    return merged


def _filter_introns_by_promoter_up(
    introns_by_chrom: Mapping[str, List[Tuple[int, int]]],
    promoter_up_by_chrom: Mapping[str, List[Tuple[int, int]]],
    max_overlaps: int,
) -> Dict[str, List[Tuple[int, int]]]:
    filtered = {}  # type: Dict[str, List[Tuple[int, int]]]
    for chrom, introns in introns_by_chrom.items():
        promoters = promoter_up_by_chrom.get(chrom, [])
        kept = []
        first_promoter = 0
        for intron_start, intron_end in introns:
            while (
                first_promoter < len(promoters)
                and promoters[first_promoter][1] <= intron_start
            ):
                first_promoter += 1
            count = 0
            index = first_promoter
            while index < len(promoters) and promoters[index][0] < intron_end:
                if promoters[index][1] > intron_start:
                    count += 1
                    if count > max_overlaps:
                        break
                index += 1
            if count <= max_overlaps:
                kept.append((intron_start, intron_end))
        if kept:
            filtered[chrom] = kept
    return filtered


def _write_merged_intervals(
    intervals_by_chrom: Mapping[str, List[Tuple[int, int]]],
    chrom_order: Tuple[str, ...],
    output_path: Path,
    merge_gap: int = 0,
) -> None:
    tmp_path = output_path.with_name(output_path.name + ".tmp")
    try:
        with tmp_path.open("w", encoding="utf-8") as handle:
            for chrom in chrom_order:
                merged = _merge_intervals(
                    list(intervals_by_chrom.get(chrom, [])),
                    max_gap=merge_gap,
                )
                for start, end in merged:
                    handle.write("{}\t{}\t{}\n".format(chrom, start, end))
        tmp_path.replace(output_path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()


def _write_intervals(
    intervals_by_chrom: Mapping[str, List[Tuple[int, int]]],
    chrom_order: Tuple[str, ...],
    output_path: Path,
) -> None:
    tmp_path = output_path.with_name(output_path.name + ".tmp")
    try:
        with tmp_path.open("w", encoding="utf-8") as handle:
            for chrom in chrom_order:
                for start, end in intervals_by_chrom.get(chrom, []):
                    handle.write("{}\t{}\t{}\n".format(chrom, start, end))
        tmp_path.replace(output_path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()


def _merge_intervals(
    intervals: List[Tuple[int, int]],
    max_gap: int = 0,
) -> List[Tuple[int, int]]:
    if not intervals:
        return []
    merged = []  # type: List[Tuple[int, int]]
    current_start, current_end = sorted(intervals)[0]
    for start, end in sorted(intervals)[1:]:
        if start - current_end <= max_gap:
            current_end = max(current_end, end)
            continue
        merged.append((current_start, current_end))
        current_start, current_end = start, end
    merged.append((current_start, current_end))
    return merged


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
    return int(float(text) * multiplier)


def _format_bp(value: int) -> str:
    if value % 1000 == 0:
        return "{}kb".format(value // 1000)
    return "{}bp".format(value)
