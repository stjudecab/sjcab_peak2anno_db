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
from ._registry import user_data_dir

PathLike = Union[str, os.PathLike]

_ATTRIBUTE_RE = re.compile(r'(\S+)\s+"([^"]*)"')
_ENSEMBL_SPECIES = {
    "human": ("homo_sapiens", "vertebrates", False, (("NCBI36", 54, 54), ("GRCh37", 55, 75), ("GRCh38", 76, 10**9))),
    "homo_sapiens": ("homo_sapiens", "vertebrates", False, (("NCBI36", 54, 54), ("GRCh37", 55, 75), ("GRCh38", 76, 10**9))),
    "mouse": ("mus_musculus", "vertebrates", False, (("NCBIM37", 54, 67), ("GRCm38", 68, 102), ("GRCm39", 103, 10**9))),
    "mus_musculus": ("mus_musculus", "vertebrates", False, (("NCBIM37", 54, 67), ("GRCm38", 68, 102), ("GRCm39", 103, 10**9))),
    "rat": ("rattus_norvegicus", "vertebrates", False, (("Rnor_6.0", 80, 104), ("mRatBN7.2", 105, 113), ("GRCr8", 114, 10**9))),
    "rattus_norvegicus": ("rattus_norvegicus", "vertebrates", False, (("Rnor_6.0", 80, 104), ("mRatBN7.2", 105, 113), ("GRCr8", 114, 10**9))),
    "dog": ("canis_familiaris", "vertebrates", False, (("CanFam3.1", 75, 10**9),)),
    "cat": ("felis_catus", "vertebrates", False, (("Felis_catus_9.0", 93, 10**9),)),
    "chicken": ("gallus_gallus", "vertebrates", False, (("Gallus_gallus-5.0", 86, 10**9),)),
    "zebrafish": ("danio_rerio", "vertebrates", False, (("GRCz11", 92, 10**9),)),
    "danio_rerio": ("danio_rerio", "vertebrates", False, (("GRCz11", 92, 10**9),)),
    "macaque": ("macaca_fascicularis", "vertebrates", False, (("Macaca_fascicularis_6.0", 103, 10**9),)),
    "rhesus": ("macaca_mulatta", "vertebrates", False, (("Mmul_10", 75, 10**9),)),
    "rabbit": ("oryctolagus_cuniculus", "vertebrates", False, (("OryCun2.0", 75, 10**9),)),
    "pig": ("sus_scrofa", "vertebrates", False, (("Sscrofa11.1", 75, 10**9),)),
    "fruitfly": ("drosophila_melanogaster", "metazoa", True, (("BDGP6.32", 103, 10**9),)),
    "drosophila": ("drosophila_melanogaster", "metazoa", True, (("BDGP6.32", 103, 10**9),)),
    "drosophila_melanogaster": ("drosophila_melanogaster", "metazoa", True, (("BDGP6.32", 103, 10**9),)),
    "worm": ("caenorhabditis_elegans", "metazoa", True, (("WBcel235", 71, 10**9),)),
    "caenorhabditis_elegans": ("caenorhabditis_elegans", "metazoa", True, (("WBcel235", 71, 10**9),)),
    "yeast": ("saccharomyces_cerevisiae", "fungi", False, (("R64-1-1", 76, 10**9),)),
    "saccharomyces_cerevisiae": ("saccharomyces_cerevisiae", "fungi", False, (("R64-1-1", 76, 10**9),)),
    "arabidopsis": ("arabidopsis_thaliana", "plants", True, (("TAIR10", 40, 10**9),)),
    "arabidopsis_thaliana": ("arabidopsis_thaliana", "plants", True, (("TAIR10", 40, 10**9),)),
    "rice": ("oryza_sativa", "plants", True, (("IRGSP-1.0", 40, 10**9),)),
    "oryza_sativa": ("oryza_sativa", "plants", True, (("IRGSP-1.0", 40, 10**9),)),
    "wheat": ("triticum_aestivum", "plants", True, (("IWGSC", 40, 10**9),)),
    "maize": ("zea_mays", "plants", True, (("Zm-B73-REFERENCE-NAM-5.0", 54, 10**9),)),
    "corn": ("zea_mays", "plants", True, (("Zm-B73-REFERENCE-NAM-5.0", 54, 10**9),)),
    "tomato": ("solanum_lycopersicum", "plants", True, (("SL3.0", 42, 10**9),)),
    "soybean": ("glycine_max", "plants", True, (("Glycine_max_v2.1", 43, 10**9),)),
    "fission_yeast": ("schizosaccharomyces_pombe", "fungi", True, (("ASM294v2", 40, 10**9),)),
    "aspergillus": ("aspergillus_nidulans", "fungi", True, (("ASM1142v1", 40, 10**9),)),
    "candida": ("candida_albicans", "fungi", True, (("GCA000182965v3", 40, 10**9),)),
    "mosquito": ("anopheles_gambiae", "metazoa", True, (("AgamP4", 40, 10**9),)),
    "plasmodium": ("plasmodium_falciparum", "protists", True, (("ASM276v2", 40, 10**9),)),
    "toxoplasma": ("toxoplasma_gondii", "protists", True, (("TGA4", 40, 10**9),)),
}
_ENSEMBL_CATALOGS = (
    (
        "vertebrates",
        False,
        "https://ftp.ebi.ac.uk/pub/ensembl/current/species_EnsemblVertebrates.txt",
        "species_EnsemblVertebrates.txt",
        "vertebrates",
    ),
    (
        "genomes",
        True,
        "https://ftp.ebi.ac.uk/pub/ensemblgenomes/current/species.txt",
        "species.txt",
        "genomes",
    ),
)
UCSC_GENES_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/{}/bigZips/genes/"


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
    metadata = _ENSEMBL_SPECIES.get(species_key)
    release = str(version).lower().replace("release-", "", 1)
    if release in {"def", "default", "current", "latest"}:
        if metadata is not None:
            release = str(_latest_ensembl_release(metadata[2]))
        else:
            metadata, release = _resolve_cached_ensembl_species(species)
    if not release.isdigit():
        raise ValueError(
            "Ensembl release must be an integer or def, got {!r}.".format(version)
        )
    release_number = int(release)
    if metadata is None:
        metadata = _resolve_cached_ensembl_species(species, release_number)[0]
    latin_name, division, genomes, references = metadata
    if genomes:
        _refresh_ensembl_default_link(release_number, genomes=True)
    reference = next(
        (name for name, first, last in references if first <= release_number <= last),
        None,
    )
    if reference is None:
        raise ValueError(
            "No Ensembl reference assembly is defined for {!r} release {}.".format(
                species, release
            )
        )
    filename = "{}.{}.{}.gtf.gz".format(latin_name.capitalize(), reference, release)
    if genomes:
        return "https://ftp.ebi.ac.uk/pub/ensemblgenomes/release-{}/{}/gtf/{}/{}".format(
            release, division, latin_name, filename
        )
    return "https://ftp.ensembl.org/pub/release-{}/gtf/{}/{}".format(
        release, latin_name, filename
    )


def ucsc_gtf_url(species: str, annotation: str = "ens") -> str:
    """Find an UCSC gene GTF for a short UCSC genome identifier."""

    listing_url = UCSC_GENES_URL.format(species)
    request = urllib.request.Request(
        listing_url, headers={"User-Agent": "sjcab_peak2anno_db"}
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        listing = response.read().decode("utf-8", "replace")
    filenames = re.findall(r'href="([^\"]+\.gtf\.gz)"', listing, re.IGNORECASE)
    if annotation.lower() in {"refseq", "ref"}:
        patterns = (".ncbiRefSeq.gtf.gz", ".refGene.gtf.gz")
    else:
        patterns = (".ensGene.gtf.gz",)
    filename = next(
        (name for name in filenames if any(name.endswith(pattern) for pattern in patterns)),
        None,
    )
    if filename is None:
        available = ", ".join(sorted(filenames)) or "none"
        raise ValueError(
            "UCSC {} annotation is unavailable for {} (available: {}).".format(
                annotation, species, available
            )
        )
    return listing_url + filename


def _resolve_cached_ensembl_species(species: str, release: Optional[int] = None):
    """Resolve an unlisted species from separate Vertebrates/Genomes catalogs."""

    wanted = species.lower().replace(" ", "_")
    root = user_data_dir() / "ensembl"
    for catalog_name, genomes, current_url, filename, cache_subdir in _ENSEMBL_CATALOGS:
        catalog_root = root / cache_subdir
        cache_path = catalog_root / filename
        if not cache_path.exists():
            cache_path = _download_ensembl_species_catalog(
                current_url, catalog_root, filename
            )
        match = _find_ensembl_species(cache_path, wanted, genomes)
        if match is None:
            continue
        if release is None:
            release = _latest_ensembl_release(genomes)
        release_path = catalog_root / str(release) / filename
        if genomes:
            release_url = "https://ftp.ebi.ac.uk/pub/ensemblgenomes/release-{}/species.txt".format(
                release
            )
            if not release_path.exists():
                release_path = _download_ensembl_species_catalog(
                    release_url, release_path.parent, filename
                )
            release_match = _find_ensembl_species(release_path, wanted, genomes)
            if release_match is not None:
                match = release_match
        else:
            release_path.parent.mkdir(parents=True, exist_ok=True)
            if not release_path.exists():
                shutil.copyfile(cache_path, release_path)
        _refresh_ensembl_species_cache_link(release_path, release, genomes)
        return match, str(release)
    raise ValueError(
        "Unknown Ensembl species {!r}; it was not found in the Ensembl catalogs.".format(
            species
        )
    )


def _download_ensembl_species_catalog(url: str, root: Path, filename: str) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "sjcab-peak2anno-db"})
    with urllib.request.urlopen(request, timeout=120) as response:
        rows = response.read().decode("utf-8", "replace").splitlines()
    destination = root / filename
    tmp_path = destination.with_name(destination.name + ".tmp")
    with tmp_path.open("w", encoding="utf-8") as handle:
        handle.write("species\tdivision\tassembly\n")
        for row in rows:
            if not row or row.startswith("#"):
                continue
            fields = row.split("\t")
            if len(fields) >= 5:
                handle.write("{}\t{}\t{}\n".format(fields[1], fields[2], fields[4]))
    tmp_path.replace(destination)
    return destination


def _find_ensembl_species(cache_path: Path, wanted: str, genomes: bool):
    with cache_path.open("r", encoding="utf-8") as handle:
        next(handle, None)
        for row in handle:
            fields = row.rstrip("\n").split("\t")
            if len(fields) == 3 and fields[0].lower() == wanted:
                division = fields[1].replace("Ensembl", "").lower()
                return fields[0], division, genomes, ((fields[2], 0, 10**9),)
    return None


def _refresh_ensembl_species_cache_link(
    cache_path: Path, release: int, genomes: bool
) -> None:
    """Expose release metadata through the conventional ``def`` cache path."""

    catalog_root = cache_path.parent.parent
    _refresh_ensembl_default_link(release, genomes=genomes)
    root_path = catalog_root / cache_path.name
    if root_path.is_symlink():
        root_path.unlink()
    if root_path.exists():
        return
    try:
        root_path.symlink_to(Path("def") / cache_path.name)
    except OSError:
        shutil.copyfile(cache_path, root_path)


def _refresh_ensembl_default_link(release: int, genomes: bool) -> None:
    """Point the matching Ensembl catalog's ``def`` link at a release."""

    ensembl_root = user_data_dir() / "ensembl"
    ensembl_root.mkdir(parents=True, exist_ok=True)
    default_dir = ensembl_root / "def"
    if default_dir.is_symlink():
        default_dir.unlink()
    elif default_dir.exists():
        try:
            default_dir.rmdir()
        except OSError:
            return
    try:
        target = Path("genomes" if genomes else "vertebrates") / str(release)
        default_dir.symlink_to(target, target_is_directory=True)
    except OSError:
        default_dir.mkdir(parents=True, exist_ok=True)


def _latest_ensembl_release(genomes: bool) -> int:
    if genomes:
        request = urllib.request.Request(
            "https://ftp.ebi.ac.uk/pub/ensemblgenomes/VERSION",
            headers={"User-Agent": "sjcab-peak2anno-db"},
        )
        with urllib.request.urlopen(request, timeout=120) as response:
            value = response.read().decode("utf-8", "replace").strip()
        match = re.search(r"\d+", value)
        if match is None:
            raise ValueError("Invalid Ensembl Genomes VERSION value: {!r}.".format(value))
        return int(match.group(0))
    server = (
        "https://ftp.ebi.ac.uk/pub/ensemblgenomes"
        if genomes
        else "https://ftp.ensembl.org"
    )
    request = urllib.request.Request(
        server + "/pub/", headers={"User-Agent": "sjcab-peak2anno-db"}
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        listing = response.read().decode("utf-8", "replace")
    releases = [int(value) for value in re.findall(r'href="release-(\d+)/', listing)]
    if not releases:
        raise ValueError("Unable to determine the latest Ensembl release.")
    return max(releases)


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
    source: str = "auto",
    cache_dir: Optional[PathLike] = None,
    ucsc_annotation: str = "ens",
) -> Path:
    """Download one GENCODE or Ensembl GTF file and return the local path."""

    output = Path(output_dir).expanduser()
    output.mkdir(parents=True, exist_ok=True)
    if gtf_url:
        url = gtf_url
    elif source.lower() == "auto":
        source = "gencode" if species.lower() in {
            "hg19", "hg38", "grch37", "grch38", "mm9", "mm10", "mm39",
        } else "ensembl"
        if source == "ensembl":
            try:
                url = ucsc_gtf_url(species, ucsc_annotation)
            except (OSError, ValueError):
                url = ensembl_gtf_url(species, version)
        else:
            url = gencode_gtf_url(species, version)
    elif source.lower() == "ensembl":
        url = ensembl_gtf_url(species, version)
    elif source.lower() == "ucsc":
        url = ucsc_gtf_url(species, ucsc_annotation)
    else:
        url = gencode_gtf_url(species, version)
    filename = url.rstrip("/").split("/")[-1]
    cache_root = user_data_dir(cache_dir) / "cachegtf"
    cache_root.mkdir(parents=True, exist_ok=True)
    destination = cache_root / filename

    existing_gtf = _find_existing_gtf(destination, output / filename)
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


def _find_existing_gtf(
    destination: Path, legacy_destination: Optional[Path] = None
) -> Optional[Path]:
    """Return an already-present GTF from output dir or current working dir."""

    if destination.exists():
        return destination
    if legacy_destination is not None and legacy_destination.exists():
        return legacy_destination

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
    source: str = "auto",
    cache_dir: Optional[PathLike] = None,
    ucsc_annotation: str = "ens",
    clean_cache: bool = False,
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
            cache_dir=cache_dir,
            ucsc_annotation=ucsc_annotation,
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

    if clean_cache and gtf_path is not None:
        cached_path = Path(gtf_path)
        cache_root = user_data_dir(cache_dir) / "cachegtf"
        if cached_path.parent == cache_root and cached_path.exists():
            cached_path.unlink()

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
