"""Download and convert GENCODE GTF files to BED resources."""

from __future__ import annotations

import gzip
import os
import re
import shutil
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Mapping, Optional, Tuple, Union

from ._derive import write_deduplong
from ._download import ProgressCallback, download_file, report_progress
from ._download_log import record_download_url

PathLike = Union[str, os.PathLike]

_ATTRIBUTE_RE = re.compile(r'(\S+)\s+"([^"]*)"')
_ENSEMBL_SPECIES = {
    "human": ("homo_sapiens", "GRCh38"),
    "homo_sapiens": ("homo_sapiens", "GRCh38"),
    "mouse": ("mus_musculus", "GRCm38"),
    "mus_musculus": ("mus_musculus", "GRCm38"),
    "rat": ("rattus_norvegicus", "Rnor_6.0"),
    "rattus_norvegicus": ("rattus_norvegicus", "Rnor_6.0"),
    "zebrafish": ("danio_rerio", "GRCz11"),
    "danio_rerio": ("danio_rerio", "GRCz11"),
    "fruitfly": ("drosophila_melanogaster", "BDGP6.32"),
    "drosophila_melanogaster": ("drosophila_melanogaster", "BDGP6.32"),
    "worm": ("caenorhabditis_elegans", "WBcel235"),
    "caenorhabditis_elegans": ("caenorhabditis_elegans", "WBcel235"),
}


@dataclass
class _TranscriptRecord:
    chrom: str
    start: int
    end: int
    gene_name: str
    strand: str
    gene_id: str
    transcript_id: str
    gene_type: str
    order: int
    exon_length: int = 0


def gencode_gtf_url(species: str, version: str) -> str:
    """Return the default public GENCODE GTF URL for a known species/version.

    ``species`` should be a genome build such as ``hg38``, ``hg19``, ``mm10``,
    or ``mm39``. Lifted human releases use versions such as ``v31lift37``.
    Older local lifted mouse builds such as ``mm9`` do not have a stable public
    GENCODE URL; pass ``gtf_url=`` or ``gtf_path=`` to
    :func:`download_and_convert_gencode_gtf` for those resources.
    """

    species_key = species.lower()
    version_key = version.lower()

    if species_key in {"hg38", "hg19", "grch38", "grch37"}:
        release = _human_release_number(version)
        base = (
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
            "release_{}/".format(release)
        )
        if "lift37" in version_key or species_key in {"hg19", "grch37"}:
            if "lift37" in version_key:
                filename_version = version
            else:
                filename_version = "{}lift37".format(version)
            return "{}GRCh37_mapping/gencode.{}.annotation.gtf.gz".format(
                base, filename_version
            )
        return "{}gencode.{}.annotation.gtf.gz".format(base, version)

    if species_key in {"mm10", "mm39", "grcm38", "grcm39"}:
        release = _mouse_release_number(version)
        return (
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/"
            "release_M{}/gencode.{}.annotation.gtf.gz".format(release, version)
        )

    raise ValueError(
        "No default GENCODE URL is available for species/build {!r}.".format(
            species
        )
    )


def ensembl_gtf_url(species: str, version: str) -> str:
    """Return an Ensembl GTF URL for a species and release."""

    species_key = species.lower().replace(" ", "_")
    latin_name, reference = _ENSEMBL_SPECIES.get(species_key, (species_key, None))
    release = str(version).lower().replace("release-", "", 1)
    if release in {"def", "default", "latest"}:
        listing_url = "https://ftp.ensembl.org/pub/current_gtf/{}/".format(latin_name)
        request = urllib.request.Request(
            listing_url, headers={"User-Agent": "sjcab-peak2anno-db"}
        )
        with urllib.request.urlopen(request, timeout=120) as response:
            listing = response.read().decode("utf-8", "replace")
        matches = re.findall(r'href="([^\"]+\\.gtf\\.gz)"', listing, re.IGNORECASE)
        if not matches:
            raise ValueError("No current Ensembl GTF found for {!r}.".format(species))
        return listing_url + matches[0]
    if not release.isdigit():
        raise ValueError(
            "Ensembl release must be an integer or def, got {!r}.".format(version)
        )
    if reference is None:
        raise ValueError(
            "Unknown Ensembl assembly for {!r}; use a known species name or --url.".format(species)
        )
    filename = "{}.{}.{}.gtf.gz".format(latin_name.capitalize(), reference, release)
    return "https://ftp.ensembl.org/pub/release-{}/gtf/{}/{}".format(
        release, latin_name, filename
    )


def gencode_bed_filename(
    species: str, version: str, include_type: bool = True
) -> str:
    """Return the conventional GENCODE BED file name for a build/version."""

    if include_type:
        suffix = "gene.bed.withtype"
    else:
        suffix = "gtf.bed"
    return "gencode.{}.{}.{}".format(version, species, suffix)


def gencode_bed_dir(output_dir: PathLike, species: str, version: str) -> Path:
    """Return the version directory for downloaded GENCODE BED resources."""

    return Path(output_dir).expanduser() / "bed" / species / version


def download_gencode_gtf(
    species: str,
    version: str,
    output_dir: PathLike = ".",
    gtf_url: Optional[str] = None,
    overwrite: bool = True,
    log_data_dir: Optional[PathLike] = None,
    progress: Optional[ProgressCallback] = None,
    source: str = "gencode",
) -> Path:
    """Download one GENCODE or Ensembl GTF file and return the local path."""

    output = Path(output_dir).expanduser()
    output.mkdir(parents=True, exist_ok=True)
    if gtf_url:
        url = gtf_url
    elif source.lower() == "ensembl":
        url = ensembl_gtf_url(species, version)
    else:
        url = gencode_gtf_url(species, version)
    filename = url.rstrip("/").split("/")[-1]
    destination = output / filename

    existing_gtf = _find_existing_gtf(destination)
    if existing_gtf is not None:
        report_progress(
            progress,
            "download-gencode-gtf: using existing {}".format(existing_gtf),
        )
        return existing_gtf

    if progress is None:
        _download_file(url, destination)
    else:
        _download_file(url, destination, progress=progress)
    record_download_url(url, data_dir=log_data_dir, destination=destination)
    return destination


def _find_existing_gtf(destination: Path) -> Optional[Path]:
    """Return an already-present GTF from output dir or current working dir."""

    if destination.exists():
        return destination

    cwd_candidate = Path.cwd() / destination.name
    if _same_path(cwd_candidate, destination):
        return None
    if cwd_candidate.exists():
        return cwd_candidate
    return None


def _same_path(left: Path, right: Path) -> bool:
    try:
        return left.resolve() == right.resolve()
    except OSError:
        return left == right


def convert_gencode_gtf_to_bed(
    gtf_path: PathLike,
    output_bed: PathLike,
    include_type: bool = True,
    gene_types: Optional[Iterable[str]] = None,
) -> Path:
    """Convert a GENCODE GTF to a transcript-level gene BED file.

    The default output matches files such as
    ``gencode.v31.hg38.gene.bed.withtype``:

    ``chrom, start, end, gene_name, transcript_exon_length, strand,``
    ``gene_id.version, transcript_id.version, gene_type``.

    Set ``include_type=False`` to write the package's compact 8-column BED
    form. ``gene_types`` can be used to keep only selected GENCODE gene types.
    """

    output = Path(output_bed).expanduser()
    output.parent.mkdir(parents=True, exist_ok=True)

    selected_types = None
    if gene_types is not None:
        selected_types = set(gene_types)

    records, order = _read_transcripts(gtf_path)
    tmp_path = output.with_name(output.name + ".tmp")
    try:
        with tmp_path.open("w", encoding="utf-8") as handle:
            for transcript_id in order:
                record = records[transcript_id]
                if selected_types is not None and record.gene_type not in selected_types:
                    continue
                length = record.exon_length or (record.end - record.start)
                fields = [
                    record.chrom,
                    str(record.start),
                    str(record.end),
                    record.gene_name,
                    str(length),
                    record.strand,
                    record.gene_id,
                    record.transcript_id,
                ]
                if include_type:
                    fields.append(record.gene_type)
                handle.write("{}\n".format("\t".join(fields)))
        tmp_path.replace(output)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()

    return output


def download_and_convert_gencode_gtf(
    species: str,
    version: str,
    output_dir: PathLike = ".",
    gtf_path: Optional[PathLike] = None,
    gtf_url: Optional[str] = None,
    output_bed: Optional[PathLike] = None,
    include_type: bool = True,
    gene_types: Optional[Iterable[str]] = None,
    overwrite: bool = True,
    log_data_dir: Optional[PathLike] = None,
    progress: Optional[ProgressCallback] = None,
    source: str = "gencode",
) -> Path:
    """Download or reuse a GENCODE GTF and convert it to BED.

    When ``output_bed`` is not provided, the generated BEDs are organized as
    ``{output_dir}/bed/{species}/{version}/{isoform_set}.{annotation}.bed``
    and ``{output_dir}/bed/{species}/def`` points at the requested version
    directory. Passing ``output_bed`` keeps the explicit single-file behavior.
    """

    target_dir = Path(output_dir).expanduser()
    target_dir.mkdir(parents=True, exist_ok=True)
    layout_mode = output_bed is None

    if gtf_path is None:
        gtf_path = download_gencode_gtf(
            species,
            version,
            target_dir,
            gtf_url=gtf_url,
            overwrite=overwrite,
            log_data_dir=log_data_dir,
            progress=progress,
            source=source,
        )
    else:
        gtf_path = Path(gtf_path).expanduser()

    if output_bed is None:
        version_dir = gencode_bed_dir(target_dir, species, version)
        output_bed = version_dir / "all.gene.bed"
    else:
        version_dir = None
        output_bed = Path(output_bed).expanduser()

    output_bed = Path(output_bed)
    if overwrite or not output_bed.exists():
        convert_gencode_gtf_to_bed(
            gtf_path,
            output_bed,
            include_type=include_type,
            gene_types=gene_types,
        )

    if layout_mode and version_dir is not None:
        _write_gencode_bed_layout(output_bed, version_dir, overwrite=overwrite)
        _refresh_default_version_dir(target_dir, species, version, version_dir)

    return output_bed


def regenerate_gencode_beds(
    output_dir: PathLike,
    specs: Iterable[Tuple[str, str, PathLike]],
    include_type: bool = True,
    overwrite: bool = True,
) -> Tuple[Path, ...]:
    """Regenerate several GENCODE BED files from local GTF paths.

    ``specs`` entries are ``(species, version, gtf_path)`` tuples. This helper
    intentionally uses local GTFs, which makes cluster regeneration reproducible
    and avoids downloading large files during package maintenance.
    """

    generated = []
    for species, version, gtf_path in specs:
        generated.append(
            download_and_convert_gencode_gtf(
                species,
                version,
                output_dir,
                gtf_path=gtf_path,
                include_type=include_type,
                overwrite=overwrite,
            )
        )
    return tuple(generated)


def _read_transcripts(
    gtf_path: PathLike,
) -> Tuple[Dict[str, _TranscriptRecord], Tuple[str, ...]]:
    records = {}  # type: Dict[str, _TranscriptRecord]
    order = []  # type: list[str]

    with _open_text(gtf_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                raise ValueError(
                    "{}:{} is not a valid GTF row".format(gtf_path, line_number)
                )
            feature = fields[2]
            if feature not in {"transcript", "exon"}:
                continue

            start = int(fields[3])
            end = int(fields[4])
            attrs = _parse_attributes(fields[8])
            transcript_id = attrs.get("transcript_id")
            gene_id = attrs.get("gene_id")
            if not transcript_id or not gene_id:
                continue

            start0 = start - 1
            gene_name = attrs.get("gene_name") or _strip_version(gene_id)
            gene_type = (
                attrs.get("gene_type")
                or attrs.get("gene_biotype")
                or attrs.get("transcript_type")
                or attrs.get("transcript_biotype")
                or "."
            )

            record = records.get(transcript_id)
            if record is None:
                record = _TranscriptRecord(
                    chrom=fields[0],
                    start=start0,
                    end=end,
                    gene_name=gene_name,
                    strand=fields[6],
                    gene_id=gene_id,
                    transcript_id=transcript_id,
                    gene_type=gene_type,
                    order=len(order),
                )
                records[transcript_id] = record
                order.append(transcript_id)
            elif feature == "transcript":
                record.chrom = fields[0]
                record.start = start0
                record.end = end
                record.gene_name = gene_name
                record.strand = fields[6]
                record.gene_id = gene_id
                record.gene_type = gene_type
            else:
                record.start = min(record.start, start0)
                record.end = max(record.end, end)

            if feature == "exon":
                record.exon_length += end - start + 1

    return records, tuple(order)


def _write_gencode_bed_layout(
    all_gene_bed: Path,
    version_dir: Path,
    overwrite: bool,
) -> None:
    deduplong_gene_bed = version_dir / "deduplong.gene.bed"
    if overwrite or not deduplong_gene_bed.exists():
        write_deduplong(all_gene_bed, deduplong_gene_bed)


def _refresh_default_version_dir(
    target_dir: Path,
    species: str,
    version: str,
    version_dir: Path,
) -> None:
    default_dir = target_dir / "bed" / species / "def"
    if default_dir.is_symlink():
        default_dir.unlink()
    elif default_dir.exists():
        if default_dir.is_dir():
            shutil.rmtree(default_dir)
        else:
            default_dir.unlink()

    try:
        default_dir.symlink_to(version, target_is_directory=True)
    except OSError:
        shutil.copytree(version_dir, default_dir, dirs_exist_ok=True)


def _parse_attributes(value: str) -> Mapping[str, str]:
    return {match.group(1): match.group(2) for match in _ATTRIBUTE_RE.finditer(value)}


def _open_text(path: PathLike):
    selected = Path(path).expanduser()
    if selected.name.endswith(".gz"):
        return gzip.open(selected, "rt", encoding="utf-8")
    return selected.open("r", encoding="utf-8")


def _download_file(
    url: str,
    destination: Path,
    progress: Optional[ProgressCallback] = None,
) -> None:
    download_file(url, destination, timeout=120, progress=progress)


def _strip_version(identifier: str) -> str:
    return identifier.split(".", 1)[0]


def _human_release_number(version: str) -> str:
    match = re.match(r"v(?P<release>\d+)", version)
    if not match:
        raise ValueError("Human GENCODE version should look like v31 or v31lift37.")
    return match.group("release")


def _mouse_release_number(version: str) -> str:
    match = re.match(r"vM(?P<release>\d+)", version)
    if not match:
        raise ValueError("Mouse GENCODE version should look like vM22.")
    return match.group("release")
