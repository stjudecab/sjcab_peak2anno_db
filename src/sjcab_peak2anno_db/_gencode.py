"""Download and convert GENCODE GTF files to BED resources."""

from __future__ import annotations

import gzip
import json
import os
import re
import shutil
import time
import urllib.parse
import urllib.request
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Mapping, Optional, Tuple, Union
from urllib.error import HTTPError

from ._derive import write_deduplong
from ._config import (
    clean_cache_files,
    configured_sizes_clean,
    configured_stale_seconds,
)
from ._download import ProgressCallback, download_file, report_progress
from ._download_log import record_download_url
from ._registry import data_root, user_data_dir

PathLike = Union[str, os.PathLike]
_ENSEMBL_CACHE_MAX_AGE = 90 * 24 * 60 * 60

_ATTRIBUTE_RE = re.compile(r'(\S+)\s+"([^"]*)"')
_ENSEMBL_SPECIES = {
    "human": ("homo_sapiens", "vertebrates", False, (("NCBI36", 54, 54), ("GRCh37", 55, 75), ("GRCh38", 76, 10**9))),
    "homo_sapiens": ("homo_sapiens", "vertebrates", False, (("NCBI36", 54, 54), ("GRCh37", 55, 75), ("GRCh38", 76, 10**9))),
    "mouse": ("mus_musculus", "vertebrates", False, (("NCBIM37", 54, 67), ("GRCm38", 68, 102), ("GRCm39", 103, 10**9))),
    "mus_musculus": ("mus_musculus", "vertebrates", False, (("NCBIM37", 54, 67), ("GRCm38", 68, 102), ("GRCm39", 103, 10**9))),
    "rat": ("rattus_norvegicus", "vertebrates", False, (("Rnor_6.0", 80, 104), ("mRatBN7.2", 105, 113), ("GRCr8", 114, 10**9))),
    "rattus_norvegicus": ("rattus_norvegicus", "vertebrates", False, (("Rnor_6.0", 80, 104), ("mRatBN7.2", 105, 113), ("GRCr8", 114, 10**9))),
    "dog": ("canis_lupus_familiaris", "vertebrates", False, (("CanFam3.1", 75, 80),("ROS_Cfam_1.0",81,10**9))),
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
_ENSEMBL_DIVISION_ORDER = {
    "vertebrates": 0,
    "plants": 1,
    "metazoa": 2,
    "fungi": 3,
    "protists": 4,
    "bacteria": 5,
}
UCSC_GENES_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/{}/bigZips/genes/"
UCSC_DOWNLOADS_URL = "https://hgdownload.soe.ucsc.edu/downloads.html"
UCSC_LIFTOVER_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/"
NCBI_EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{}"


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


def ensembl_gtf_url(
    species: str,
    version: str,
    cache_dir: Optional[PathLike] = None,
    use_catalog: bool = False,
) -> str:
    """Return an Ensembl GTF URL for a species and release."""

    species_key = species.lower().replace(" ", "_")
    metadata = None if use_catalog else _ENSEMBL_SPECIES.get(species_key)
    release = str(version).lower().replace("release-", "", 1)
    if release in {"def", "default", "current", "latest"}:
        if metadata is not None:
            release = str(_latest_ensembl_release(metadata[2], cache_dir=cache_dir))
        else:
            metadata, release = _resolve_cached_ensembl_species(
                species, cache_dir=cache_dir
            )
    if not release.isdigit():
        raise ValueError(
            "Ensembl release must be an integer or def, got {!r}.".format(version)
        )
    release_number = int(release)
    if metadata is None:
        metadata = _resolve_cached_ensembl_species(
            species, release_number, cache_dir=cache_dir
        )[0]
    latin_name, division, genomes, references = metadata
    _refresh_ensembl_default_link(
        release_number, genomes=genomes, cache_dir=cache_dir
    )
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


def ensembl_assembly_name(
    species: str, version: str, cache_dir: Optional[PathLike] = None
) -> str:
    """Return the Ensembl assembly name used for a species and release."""

    species_key = species.lower().replace(" ", "_")
    metadata = _ENSEMBL_SPECIES.get(species_key)
    release = str(version).lower().replace("release-", "", 1)
    if release in {"def", "default", "current", "latest"}:
        if metadata is not None:
            release = str(_latest_ensembl_release(metadata[2], cache_dir=cache_dir))
        else:
            metadata, release = _resolve_cached_ensembl_species(
                species, cache_dir=cache_dir
            )
    if metadata is None:
        metadata = _resolve_cached_ensembl_species(
            species, int(release), cache_dir=cache_dir
        )[0]
    release_number = int(release)
    return next(
        (name for name, first, last in metadata[3] if first <= release_number <= last),
        species,
    )


def data_species_name(
    species: str, version: str, cache_dir: Optional[PathLike] = None
) -> str:
    """Return the storage name, using an assembly only for Ensembl data."""

    if species.lower() in {
        "hg19", "hg38", "grch37", "grch38", "mm9", "mm10", "mm39"
    }:
        return species
    if _match_ucsc_build(species, _ucsc_gtf_builds(cache_dir)) is not None:
        return _match_ucsc_build(species, _ucsc_gtf_builds(cache_dir))
    return ensembl_assembly_name(species, version, cache_dir=cache_dir)


def ensembl_assembly_accession(
    species: str, version: str, cache_dir: Optional[PathLike] = None
) -> Optional[str]:
    """Return a cached Ensembl assembly accession when one is available."""

    root = user_data_dir(cache_dir) / "ensembl"
    species_key = species.lower().replace(" ", "_")
    search_species = _ENSEMBL_SPECIES.get(species_key, (species,))[0]
    for _catalog_name, _genomes, _url, filename, cache_subdir in _ENSEMBL_CATALOGS:
        path = root / cache_subdir / "def" / filename
        legacy_path = root / cache_subdir / filename
        if not path.exists():
            try:
                _resolve_cached_ensembl_species(species, cache_dir=cache_dir)
            except (OSError, ValueError):
                pass
        if not path.exists() and legacy_path.exists():
            path = legacy_path
        if not path.exists():
            continue
        rows = _read_ensembl_species_catalog(path)
        rows.sort(key=lambda fields: _ensembl_division_rank(fields[2]))
        wanted = _normalize_species_match(search_species)
        for fields in rows:
            if any(
                _normalize_species_match(fields[index]) == wanted
                for index in (0, 1, 3)
            ):
                return fields[4] or None
    return None


def ucsc_gtf_url(
    species: str,
    annotation: str = "ens",
    cache_dir: Optional[PathLike] = None,
) -> str:
    """Find an UCSC gene GTF for a short UCSC genome identifier."""

    build = _match_ucsc_build(species, _ucsc_gtf_builds(cache_dir))
    if build is None:
        raise ValueError(
            "UCSC build {!r} is not listed with a GTF on {}.".format(
                species, UCSC_DOWNLOADS_URL
            )
        )
    listing_url = UCSC_GENES_URL.format(build)
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


def _match_ucsc_build(species: str, builds) -> Optional[str]:
    """Return the canonical UCSC build ID for a normalized request."""

    wanted = _normalize_species_match(species)
    return next(
        (build for build in builds if _normalize_species_match(build) == wanted),
        None,
    )


def _ucsc_gtf_builds(
    cache_dir: Optional[PathLike] = None,
    sizes_clean: Optional[bool] = None,
):
    """Return UCSC build IDs, refreshing the runtime cache when stale."""

    cache_path = user_data_dir(cache_dir) / "ucsc" / "gtf_builds.tsv"
    old_builds = _read_ucsc_builds(cache_path) if cache_path.exists() else set()
    if cache_path.exists() and time.time() - cache_path.stat().st_mtime < configured_stale_seconds():
        if _ucsc_chain_rows_need_migration(cache_path) or not _read_ucsc_chain_urls(cache_path):
            _cache_ucsc_chain_urls(cache_path)
        return old_builds
    package_path = data_root() / "gtf_builds.tsv"
    if not cache_path.exists() and package_path.exists():
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(package_path, cache_path)
        _cache_ucsc_chain_urls(cache_path)
        return _read_ucsc_builds(cache_path)

    try:
        request = urllib.request.Request(
            UCSC_DOWNLOADS_URL, headers={"User-Agent": "sjcab_peak2anno_db"}
        )
        with urllib.request.urlopen(request, timeout=120) as response:
            page = response.read().decode("utf-8", "replace")
    except OSError:
        if cache_path.exists():
            return _read_ucsc_builds(cache_path)
        if package_path.exists():
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(package_path, cache_path)
            _cache_ucsc_chain_urls(cache_path)
            return _read_ucsc_builds(cache_path)
        raise

    builds = set()
    metadata_path = cache_path if cache_path.exists() else package_path
    descriptions = _read_ucsc_descriptions(metadata_path) if metadata_path.exists() else {}
    hrefs = re.findall(r'href=["\']([^"\']+)["\']', page, re.IGNORECASE)
    for href in hrefs:
        match = re.search(
            r"(?:^|/)goldenPath/([A-Za-z0-9_.-]+)(?:/[^?#]*)(?:genes|\.gtf(?:\.gz)?)",
            href,
            re.IGNORECASE,
        )
        if match is None:
            match = re.search(
                r"(?:^|/)goldenPath/([A-Za-z0-9_.-]+)/bigZips/",
                href,
                re.IGNORECASE,
            )
        if match:
            builds.add(match.group(1))

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = cache_path.with_name(cache_path.name + ".tmp")
    chain_urls = _download_ucsc_chain_urls()
    with tmp_path.open("w", encoding="utf-8") as handle:
        for build in sorted(builds):
            description = descriptions.get(build)
            handle.write(
                "{}\t{}\n".format(build, description)
                if description
                else build + "\n"
            )
    _merge_ucsc_chain_rows(tmp_path, chain_urls)
    tmp_path.replace(cache_path)
    _cache_new_ucsc_sizes(builds - old_builds, cache_dir, sizes_clean=sizes_clean)
    return builds


def _read_ucsc_builds(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        return {
            row.split("\t", 1)[0].strip()
            for row in handle
            if row.strip() and not row.startswith("#")
    }


def _read_ucsc_descriptions(path: Path):
    descriptions = {}
    with path.open("r", encoding="utf-8") as handle:
        for row in handle:
            fields = row.rstrip("\n").split("\t")
            if len(fields) >= 2 and fields[0].strip() and not fields[0].startswith("#"):
                descriptions[fields[0].strip()] = fields[1].strip()
    return descriptions


def _read_ucsc_chain_urls(path: Path):
    urls = {}
    with path.open("r", encoding="utf-8") as handle:
        for row in handle:
            fields = row.rstrip("\n").split("\t")
            if len(fields) >= 4 and fields[2].strip().lower().endswith(".over.chain.gz"):
                urls[fields[2].strip()] = fields[3].strip()
            elif len(fields) >= 3 and fields[0].strip().lower() == "#liftover":
                urls[fields[1].strip()] = fields[2].strip()
    return urls


def _ucsc_chain_rows_need_migration(path: Path) -> bool:
    with path.open("r", encoding="utf-8") as handle:
        return any(row.lower().startswith("#liftover\t") for row in handle)


def _download_ucsc_chain_urls():
    try:
        request = urllib.request.Request(
            UCSC_LIFTOVER_URL, headers={"User-Agent": "sjcab_peak2anno_db"}
        )
        with urllib.request.urlopen(request, timeout=120) as response:
            page = response.read().decode("utf-8", "replace")
    except OSError:
        return {}
    urls = {}
    for href in re.findall(r'href=["\']([^"\']+\.over\.chain\.gz)["\']', page, re.IGNORECASE):
        url = urllib.parse.urljoin(UCSC_LIFTOVER_URL, href)
        urls[Path(urllib.parse.urlparse(url).path).name] = url
    return urls


def _cache_ucsc_chain_urls(cache_path: Path):
    urls = _download_ucsc_chain_urls()
    if not urls:
        return
    _merge_ucsc_chain_rows(cache_path, urls)


def _merge_ucsc_chain_rows(path: Path, chain_urls: Dict[str, str]) -> None:
    """Attach cached liftOver chain names and URLs to UCSC build rows."""

    rows = []
    for row in path.read_text(encoding="utf-8").splitlines():
        fields = row.split("\t")
        if not row.strip() or row.startswith("#"):
            continue
        build = fields[0].strip()
        if not build:
            continue
        chain = _chain_for_ucsc_build(build, chain_urls)
        if chain is not None:
            fields = fields[:2] + [chain[0], chain[1]]
        rows.append("\t".join(fields))
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def _chain_for_ucsc_build(build: str, chain_urls: Dict[str, str]):
    wanted = build.casefold()
    matches = []
    for name, url in chain_urls.items():
        match = re.fullmatch(r"hg38to(.+?)\.over\.chain\.gz", name, re.IGNORECASE)
        if match is None:
            continue
        target = match.group(1)
        if target.casefold() == wanted:
            return name, url
        target_base = re.sub(r"[0-9]+$", "", target).casefold()
        build_base = re.sub(r"[0-9]+$", "", build).casefold()
        if target_base == build_base:
            suffix = re.search(r"[0-9]+$", target)
            matches.append((int(suffix.group()) if suffix else -1, name, url))
    if matches:
        _, name, url = max(matches)
        return name, url
    return None


def ucsc_liftover_chain_url(
    source_genome: str,
    chain_name: str,
    cache_dir: Optional[PathLike] = None,
) -> str:
    """Resolve a hg38 liftOver chain from the cached UCSC directory listing."""

    if str(source_genome).lower() == "hg38":
        cache_path = user_data_dir(cache_dir) / "ucsc" / "gtf_builds.tsv"
        _ucsc_gtf_builds(cache_dir)
        requested = re.fullmatch(
            r"(.+?)to(.+?)\.over\.chain\.gz", chain_name, re.IGNORECASE
        )
        for name, url in _read_ucsc_chain_urls(cache_path).items():
            cached = re.fullmatch(
                r"(.+?)to(.+?)\.over\.chain\.gz", name, re.IGNORECASE
            )
            if cached and requested and tuple(
                part.casefold() for part in cached.groups()
            ) == tuple(part.casefold() for part in requested.groups()):
                return url
    return urllib.parse.urljoin(
        "https://hgdownload.soe.ucsc.edu/goldenPath/{}/liftOver/".format(source_genome),
        chain_name,
    )


def _cache_new_ucsc_sizes(
    builds, cache_dir: Optional[PathLike], sizes_clean: Optional[bool] = None
) -> None:
    sizes_dir = user_data_dir(cache_dir) / "sizes"
    sizes_dir.mkdir(parents=True, exist_ok=True)
    for build in sorted(builds):
        raw_path = sizes_dir / "{}.sizes".format(build)
        clean_path = sizes_dir / "{}.sizes.clean".format(build)
        effective_sizes_clean = (
            configured_sizes_clean() if sizes_clean is None else sizes_clean
        )
        if raw_path.exists() and (not effective_sizes_clean or clean_path.exists()):
            continue
        url = "https://hgdownload.soe.ucsc.edu/goldenPath/{0}/bigZips/{0}.chrom.sizes".format(
            build
        )
        if not effective_sizes_clean:
            continue
        try:
            request = urllib.request.Request(
                url, headers={"User-Agent": "sjcab_peak2anno_db"}
            )
            with urllib.request.urlopen(request, timeout=120) as response:
                rows = response.read().decode("utf-8").splitlines()
        except OSError:
            continue
        temporary = raw_path.with_name(raw_path.name + ".tmp")
        with temporary.open("w", encoding="utf-8") as handle:
            for row in rows:
                fields = row.split()
                if len(fields) >= 2:
                    handle.write("{}\t{}\n".format(fields[0], fields[1]))
        temporary.replace(raw_path)
        try:
            primary_chromosomes = _download_ucsc_primary_chromosomes(build)
        except (OSError, ValueError, KeyError, json.JSONDecodeError):
            continue
        temporary = clean_path.with_name(clean_path.name + ".tmp")
        with temporary.open("w", encoding="utf-8") as handle:
            for row in rows:
                fields = row.split()
                if len(fields) >= 2 and fields[0] in primary_chromosomes:
                    handle.write("{}\t{}\n".format(fields[0], fields[1]))
        temporary.replace(clean_path)


def _download_ucsc_primary_chromosomes(species: str):
    """Return assembled molecules from the NCBI report for a UCSC build."""

    report_cache = user_data_dir() / "cache" / "assembly_reports" / (
        species + ".assembly_report.txt"
    )
    if report_cache.exists():
        return _parse_assembly_report(
            report_cache.read_text(encoding="utf-8"), species
        )

    find_url = "https://api.genome.ucsc.edu/findGenome?q={}".format(
        urllib.parse.quote(species)
    )
    find_payload = json.loads(
        urllib.request.urlopen(find_url, timeout=120).read().decode("utf-8")
    )
    genome = find_payload.get(species)
    if genome is None:
        genome = next(
            (value for key, value in find_payload.items() if key.lower() == species.lower()),
            None,
        )
    if not genome or not genome.get("scientificName"):
        raise ValueError("UCSC findGenome found no scientific name for {!r}".format(species))
    record_download_url(find_url, destination=report_cache)
    scientific_name = genome["scientificName"]
    accessions = re.findall(r"(?:GC[AF]_\d+\.\d+)", genome.get("description", ""))

    ids = []
    search_terms = [
        '{} AND "{}"[Organism]'.format(species, scientific_name),
        "{}[Assembly Name]".format(species),
    ]
    search_terms.extend("{}[Assembly Accession]".format(accession) for accession in accessions)
    for search_term in search_terms:
        search_url = "{}?db=assembly&term={}".format(
            NCBI_EUTILS_URL.format("esearch.fcgi"),
            urllib.parse.quote(search_term),
        )
        search_payload = urllib.request.urlopen(search_url, timeout=120).read().decode(
            "utf-8", "replace"
        )
        record_download_url(search_url, destination=report_cache)
        ids = re.findall(r"<Id>(\d+)</Id>", search_payload)
        if ids:
            break
    if not ids:
        raise ValueError("NCBI Assembly search found no record for {!r}".format(species))
    summary_url = "{}?db=assembly&id={}&retmode=json".format(
        NCBI_EUTILS_URL.format("esummary.fcgi"), ",".join(ids[:20])
    )
    summary = json.loads(
        urllib.request.urlopen(summary_url, timeout=120).read().decode("utf-8")
    )
    record_download_url(summary_url, destination=report_cache)
    documents = summary.get("result", {})
    document = documents.get(ids[0], {})
    for candidate_id in summary.get("result", {}).get("uids", []):
        candidate = documents.get(candidate_id, {})
        if candidate.get("ucscname", "").lower() == species.lower():
            document = candidate
            break
    ftp_path = document.get("ftppath_assembly_rpt")
    if ftp_path:
        report_url = ftp_path.replace("ftp://", "https://")
    else:
        ftp_path = document.get("ftppath_refseq") or document.get("ftppath_genbank")
    if not ftp_path and document.get("assemblyaccession"):
        accession = document["assemblyaccession"]
        accession_number = accession.split("_", 1)[1].split(".", 1)[0]
        assembly_name = document.get("assemblyname", accession)
        ftp_path = "https://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}/{}/{}/{}_{}".format(
            accession[:3],
            accession_number[:3],
            accession_number[3:6],
            accession_number[6:9],
            accession,
            assembly_name,
        )
    if not ftp_path:
        raise ValueError("NCBI Assembly summary has no FTP path for {!r}".format(species))
    if not document.get("ftppath_assembly_rpt"):
        ftp_path = ftp_path.replace("ftp://", "https://")
        report_url = ftp_path.rstrip("/") + "/" + ftp_path.rstrip("/").rsplit("/", 1)[-1]
        report_url += "_assembly_report.txt"
    report = urllib.request.urlopen(report_url, timeout=120).read().decode(
        "utf-8", "replace"
    )
    report_cache.parent.mkdir(parents=True, exist_ok=True)
    report_cache.write_text(report, encoding="utf-8")
    record_download_url(report_url, destination=report_cache)
    return _parse_assembly_report(report, species)


def _parse_assembly_report(report: str, species: str):
    primary = set()
    for line in report.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) > 1 and fields[1] == "assembled-molecule":
            primary.update(
                field
                for index in (0, 2, 9)
                if len(fields) > index
                for field in (fields[index],)
                if field and field != "na"
            )
    if not primary:
        raise ValueError("NCBI Assembly report has no assembled molecules for {!r}".format(species))
    return primary


def _resolve_cached_ensembl_species(
    species: str,
    release: Optional[int] = None,
    cache_dir: Optional[PathLike] = None,
):
    """Resolve an unlisted species from separate Vertebrates/Genomes catalogs."""

    species_key = species.lower().replace(" ", "_")
    wanted = _ENSEMBL_SPECIES.get(species_key, (species_key,))[0]
    root = user_data_dir(cache_dir) / "ensembl"
    requested_release = release
    for catalog_name, genomes, current_url, filename, cache_subdir in _ENSEMBL_CATALOGS:
        catalog_root = root / cache_subdir
        catalog_release = requested_release
        if catalog_release is None:
            catalog_release = _latest_ensembl_release(genomes, cache_dir=cache_dir)
        release_path = catalog_root / str(catalog_release) / filename
        cache_path = release_path
        if not cache_path.exists() or _is_stale(cache_path):
            legacy_path = catalog_root / filename
            if legacy_path.exists():
                release_path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(legacy_path, release_path)
                cache_path = release_path
            if cache_path.exists() and not _is_stale(cache_path):
                pass
            else:
                cache_path = _download_ensembl_species_catalog(
                    _ensembl_release_catalog_url(
                        catalog_release, genomes
                    ),
                    release_path.parent,
                    filename,
                    cache_dir=cache_dir,
                )
        match = _find_ensembl_species(cache_path, wanted, genomes)
        if match is None:
            continue
        _refresh_ensembl_species_cache_link(
            cache_path, catalog_release, genomes, cache_dir=cache_dir
        )
        return match, str(catalog_release)
    candidates = []
    for _catalog_name, _genomes, _current_url, filename, cache_subdir in _ENSEMBL_CATALOGS:
        catalog_root = root / cache_subdir
        catalog_paths = [catalog_root / filename]
        catalog_paths.extend(catalog_root.glob("*/" + filename))
        for catalog_path in catalog_paths:
            if catalog_path.exists():
                candidates.extend(_find_ensembl_species_candidates(catalog_path, species))
    message = [
        "Unknown Ensembl species {!r}; it was not found in the catalogs.".format(
            species
        )
    ]
    if candidates:
        message.append("Candidates (copy column 1 as the species ID):")
        message.append("assembly\tspecies\tdivision\tname\tassembly_accession")
        message.extend(
            "\t".join((row[3], row[0], row[2], row[1], row[4]))
            for row in candidates[:20]
        )
    message.append(
        "Reminder: older Ensembl releases were not searched. Try "
        "https://www.ncbi.nlm.nih.gov/datasets/genome/ for older annotations."
    )
    raise ValueError("\n".join(message))


def _ensembl_release_catalog_url(
    release: int, genomes: bool
) -> str:
    if genomes:
        return "https://ftp.ebi.ac.uk/pub/ensemblgenomes/release-{}/species.txt".format(
            release
        )
    return "https://ftp.ensembl.org/pub/release-{}/species_EnsemblVertebrates.txt".format(
        release
    )


def _download_ensembl_species_catalog(
    url: str,
    root: Path,
    filename: str,
    cache_dir: Optional[PathLike] = None,
) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "sjcab-peak2anno-db"})
    with urllib.request.urlopen(request, timeout=120) as response:
        rows = response.read().decode("utf-8", "replace").splitlines()
    destination = root / filename
    tmp_path = destination.with_name(destination.name + ".tmp")
    with tmp_path.open("w", encoding="utf-8") as handle:
        handle.write("assembly\tspecies\tdivision\tname\tassembly_accession\n")
        for row in rows:
            if not row or row.startswith("#"):
                continue
            fields = row.split("\t")
            if len(fields) >= 5:
                fields = [field.replace(" ", "_") for field in fields]
                handle.write(
                    "{}\t{}\t{}\t{}\t{}\n".format(
                        fields[4],
                        fields[1],
                        fields[2],
                        fields[0],
                        fields[5] if len(fields) > 5 else "",
                    )
                )
    tmp_path.replace(destination)
    record_download_url(url, data_dir=cache_dir, destination=destination)
    return destination


def _find_ensembl_species(cache_path: Path, wanted: str, genomes: bool):
    wanted_normalized = _normalize_species_match(wanted)
    rows = _read_ensembl_species_catalog(cache_path)
    rows.sort(key=lambda fields: _ensembl_division_rank(fields[2]))
    for field_index in (3, 0, 1):
        for fields in rows:
            if _normalize_species_match(fields[field_index]) == wanted_normalized:
                return _ensembl_species_metadata(fields, genomes)
    return None


def _read_ensembl_species_catalog(cache_path: Path):
    rows = []
    with cache_path.open("r", encoding="utf-8") as handle:
        header = next(handle, "")
        header_fields = header.rstrip("\n").split("\t")
        new_format = len(header_fields) >= 5 and header_fields[0].lower() == "assembly"
        for row in handle:
            fields = row.rstrip("\n").split("\t")
            if new_format and len(fields) >= 5:
                rows.append((fields[1], fields[3], fields[2], fields[0], fields[4]))
            elif len(fields) >= 6:
                # Original Ensembl catalog: name, species, division, taxon,
                # assembly, assembly_accession, ...
                rows.append((fields[1], fields[0], fields[2], fields[4], fields[5]))
            elif len(fields) == 4:
                # Cache format written before assembly_accession was added.
                if fields[0].lower() == "assembly":
                    rows.append((fields[1], fields[3], fields[2], "", ""))
                else:
                    rows.append((fields[0], fields[1], fields[2], fields[3], ""))
            elif len(fields) == 3:
                rows.append((fields[0], "", fields[1], fields[2], ""))
    return rows


def _ensembl_species_metadata(fields, genomes: bool):
    species, _name, division, assembly = fields[:4]
    latin_name = species or _name
    return latin_name, division.replace("Ensembl", "").lower(), genomes, (
        (assembly, 0, 10**9),
    )


def _is_http_404(error: BaseException) -> bool:
    current = error
    while current is not None:
        if isinstance(current, HTTPError) and current.code == 404:
            return True
        current = current.__cause__
    return False


def _find_ensembl_species_candidates(cache_path: Path, wanted: str):
    wanted_normalized = _normalize_species_match(wanted)
    if not wanted_normalized:
        return []
    candidates = []
    for fields in _read_ensembl_species_catalog(cache_path):
        values = tuple(_normalize_species_match(value) for value in fields)
        if any(
            value
            and (wanted_normalized in value or value in wanted_normalized)
            for value in values
        ):
            candidates.append(fields)
    candidates.sort(key=lambda fields: _ensembl_division_rank(fields[2]))
    return candidates


def _normalize_species_match(value: str) -> str:
    return " ".join(re.sub(r"[^a-z0-9]+", " ", str(value).lower()).split())


def _ensembl_division_rank(division: str) -> int:
    normalized = _normalize_species_match(division)
    normalized = normalized.removeprefix("ensembl ")
    return _ENSEMBL_DIVISION_ORDER.get(normalized, len(_ENSEMBL_DIVISION_ORDER))


def _refresh_ensembl_species_cache_link(
    cache_path: Path,
    release: int,
    genomes: bool,
    cache_dir: Optional[PathLike] = None,
) -> None:
    """Expose release metadata through the conventional ``def`` cache path."""

    _refresh_ensembl_default_link(release, genomes=genomes, cache_dir=cache_dir)


def _refresh_ensembl_default_link(
    release: int, genomes: bool, cache_dir: Optional[PathLike] = None
) -> None:
    """Point the matching Ensembl catalog's ``def`` link at a release."""

    catalog_root = user_data_dir(cache_dir) / "ensembl" / (
        "genomes" if genomes else "vertebrates"
    )
    catalog_root.mkdir(parents=True, exist_ok=True)
    (catalog_root / str(release)).mkdir(parents=True, exist_ok=True)
    default_dir = catalog_root / "def"
    if default_dir.is_symlink():
        default_dir.unlink()
    elif default_dir.exists():
        try:
            default_dir.rmdir()
        except OSError:
            release_dir = catalog_root / str(release)
            for child in default_dir.iterdir():
                target = release_dir / child.name
                if not target.exists():
                    shutil.move(str(child), str(target))
            default_dir.rmdir()
    try:
        target = Path(str(release))
        default_dir.symlink_to(target, target_is_directory=True)
    except OSError:
        default_dir.mkdir(parents=True, exist_ok=True)


def _latest_ensembl_release(
    genomes: bool, cache_dir: Optional[PathLike] = None
) -> int:
    root = user_data_dir(cache_dir) / "ensembl" / (
        "genomes" if genomes else "vertebrates"
    )
    version_path = root / "VERSION"
    if version_path.is_file() and not _is_stale(version_path):
        value = version_path.read_text(encoding="utf-8").strip()
    else:
        version_url = (
            "https://ftp.ebi.ac.uk/pub/ensemblgenomes/VERSION"
            if genomes
            else "https://ftp.ebi.ac.uk/pub/ensembl/VERSION"
        )
        request = urllib.request.Request(
            version_url, headers={"User-Agent": "sjcab-peak2anno-db"}
        )
        with urllib.request.urlopen(request, timeout=120) as response:
            value = response.read().decode("utf-8", "replace").strip()
        root.mkdir(parents=True, exist_ok=True)
        temporary = version_path.with_name(version_path.name + ".tmp")
        temporary.write_text(value + "\n", encoding="utf-8")
        temporary.replace(version_path)
    match = re.search(r"\d+", value)
    if match is None:
        source = "Ensembl Genomes" if genomes else "Ensembl"
        raise ValueError("Invalid {} VERSION value: {!r}.".format(source, value))
    return int(match.group(0))


def _is_stale(path: Path) -> bool:
    return time.time() - path.stat().st_mtime > configured_stale_seconds()


def gencode_bed_filename(species: str, version: str) -> str:
    """Return the conventional GENCODE BED file name for a build/version."""

    return "gencode.{}.{}.gene.bed.withtype".format(version, species)


def gencode_bed_dir(output_dir: PathLike, species: str, version: str) -> Path:
    """Return the version directory for downloaded GENCODE BED resources."""

    return Path(output_dir).expanduser() / "bed" / species / version


def resolve_gencode_gtf_url(
    species: str,
    version: str,
    gtf_url: Optional[str] = None,
    source: str = "auto",
    cache_dir: Optional[PathLike] = None,
    ucsc_annotation: str = "ens",
) -> str:
    """Resolve the downloadable GTF URL without downloading it."""

    if gtf_url:
        return gtf_url
    source_key = source.lower()
    if source_key == "auto":
        source_key = (
            "gencode"
            if species.lower()
            in {"hg19", "hg38", "grch37", "grch38", "mm9", "mm10", "mm39"}
            else "ensembl"
        )
    if source_key == "ensembl":
        try:
            return ucsc_gtf_url(species, ucsc_annotation, cache_dir=cache_dir)
        except (OSError, ValueError):
            return ensembl_gtf_url(species, version, cache_dir=cache_dir)
    if source_key == "ucsc":
        return ucsc_gtf_url(species, ucsc_annotation, cache_dir=cache_dir)
    return gencode_gtf_url(species, version)


def gtf_url_is_available(url: str) -> bool:
    """Return whether a downloadable GTF URL responds without HTTP 404."""

    request = urllib.request.Request(url, headers={"User-Agent": "sjcab-peak2anno-db"})
    try:
        response = urllib.request.urlopen(request, timeout=120)
        response.close()
        return True
    except HTTPError as error:
        if error.code == 404:
            return False
        raise


def ensembl_gtf_url_candidates(
    species: str, version: str, cache_dir: Optional[PathLike] = None
):
    """Return exact catalog matches and their candidate GTF URLs."""

    release_text = str(version).lower().replace("release-", "", 1)
    wanted = _ENSEMBL_SPECIES.get(
        species.lower().replace(" ", "_"), (species,)
    )[0]
    candidates = []
    seen = set()
    for _catalog_name, genomes, _current_url, filename, cache_subdir in _ENSEMBL_CATALOGS:
        release = (
            _latest_ensembl_release(genomes, cache_dir=cache_dir)
            if release_text in {"def", "default", "current", "latest"}
            else int(release_text)
        )
        catalog_root = user_data_dir(cache_dir) / "ensembl" / cache_subdir
        catalog_path = catalog_root / str(release) / filename
        if not catalog_path.exists() or _is_stale(catalog_path):
            catalog_path = _download_ensembl_species_catalog(
                _ensembl_release_catalog_url(release, genomes),
                catalog_path.parent,
                filename,
                cache_dir=cache_dir,
            )
        rows = _read_ensembl_species_catalog(catalog_path)
        rows.sort(key=lambda fields: _ensembl_division_rank(fields[2]))
        for field_index in (3, 0, 1):
            matches = [
                fields
                for fields in rows
                if _normalize_species_match(fields[field_index])
                == _normalize_species_match(wanted)
            ]
            if not matches:
                continue
            for fields in matches:
                key = (fields[0], fields[1], fields[2], fields[3], fields[4])
                if key in seen:
                    continue
                seen.add(key)
                metadata = _ensembl_species_metadata(fields, genomes)
                latin_name, division, is_genomes, references = metadata
                reference = next(
                    (
                        name
                        for name, first, last in references
                        if first <= release <= last
                    ),
                    fields[3],
                )
                filename_url = "{}.{}.{}.gtf.gz".format(
                    latin_name.capitalize(), reference, release
                )
                if is_genomes:
                    url = "https://ftp.ebi.ac.uk/pub/ensemblgenomes/release-{}/{}/gtf/{}/{}".format(
                        release, division, latin_name, filename_url
                    )
                else:
                    url = "https://ftp.ensembl.org/pub/release-{}/gtf/{}/{}".format(
                        release, latin_name, filename_url
                    )
                candidates.append((fields, url))
            break
    return candidates


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
    selected_source = source.lower()
    if selected_source == "auto" and species.lower() not in {
        "hg19", "hg38", "grch37", "grch38", "mm9", "mm10", "mm39",
    }:
        selected_source = "ensembl"
    url = resolve_gencode_gtf_url(
        species,
        version,
        gtf_url=gtf_url,
        source=source,
        cache_dir=cache_dir,
        ucsc_annotation=ucsc_annotation,
    )
    filename = url.rstrip("/").split("/")[-1]
    cache_root = user_data_dir(cache_dir) / "cache"
    cache_root.mkdir(parents=True, exist_ok=True)
    destination = cache_root / filename

    existing_gtf = _find_existing_gtf(destination, output / filename)
    if existing_gtf is not None:
        report_progress(
            progress,
            "download-gtf: using existing {}".format(existing_gtf),
        )
        return existing_gtf

    try:
        if progress is None:
            _download_file(url, destination)
        else:
            _download_file(url, destination, progress=progress)
    except RuntimeError as error:
        if not (
            selected_source in {"auto", "ensembl"}
            and _is_http_404(error)
        ):
            raise
        report_progress(
            progress,
            "download-gtf: Ensembl URL returned 404; retrying via species catalog",
        )
        url = ensembl_gtf_url(
            species,
            version,
            cache_dir=cache_dir,
            use_catalog=True,
        )
        filename = url.rstrip("/").split("/")[-1]
        destination = cache_root / filename
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
    gene_types: Optional[Iterable[str]] = None,
    processes: int = 1,
) -> Path:
    """Convert a GENCODE GTF to a transcript-level gene BED file.

    The default output matches files such as
    ``gencode.v31.hg38.gene.bed.withtype``:

    ``chrom, start, end, gene_name, transcript_exon_length, strand,``
    ``gene_id.version, transcript_id.version, gene_type``.

    ``gene_types`` can be used to keep only selected GENCODE gene types.
    ``processes`` controls parallel transcript formatting; ``1`` keeps the
    single-process behavior.
    """

    output = Path(output_bed).expanduser()
    output.parent.mkdir(parents=True, exist_ok=True)

    if processes < 1:
        raise ValueError("processes must be at least 1")
    selected_types = None
    if gene_types is not None:
        selected_types = set(gene_types)

    records, order = _read_transcripts(gtf_path)
    tmp_path = output.with_name(output.name + ".tmp")
    try:
        records_in_order = [records[transcript_id] for transcript_id in order]
        if processes > 1:
            with ProcessPoolExecutor(max_workers=processes) as executor:
                lines = executor.map(
                    _format_transcript_bed,
                    records_in_order,
                    [selected_types] * len(records_in_order),
                    chunksize=max(1, len(records_in_order) // (processes * 4) or 1),
                )
                formatted_lines = lines
        else:
            formatted_lines = (
                _format_transcript_bed(record, selected_types)
                for record in records_in_order
            )
        with tmp_path.open("w", encoding="utf-8") as handle:
            for line in formatted_lines:
                if line:
                    handle.write(line)
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
    gene_types: Optional[Iterable[str]] = None,
    overwrite: bool = True,
    log_data_dir: Optional[PathLike] = None,
    progress: Optional[ProgressCallback] = None,
    source: str = "auto",
    cache_dir: Optional[PathLike] = None,
    ucsc_annotation: str = "ens",
    clean_cache: Optional[Union[bool, int]] = None,
    processes: int = 1,
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
    storage_species = species
    if gtf_path is None and source.lower() in {"auto", "ensembl"}:
        storage_species = data_species_name(species, version, cache_dir=cache_dir)

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
        version_dir = gencode_bed_dir(target_dir, storage_species, version)
        output_bed = version_dir / "all.gene.bed"
    else:
        version_dir = None
        output_bed = Path(output_bed).expanduser()

    output_bed = Path(output_bed)
    if overwrite or not output_bed.exists():
        convert_gencode_gtf_to_bed(
            gtf_path,
            output_bed,
            gene_types=gene_types,
            processes=processes,
        )

    if layout_mode and version_dir is not None:
        _write_gencode_bed_layout(output_bed, version_dir, overwrite=overwrite)
        _refresh_default_version_dir(target_dir, species, version, version_dir)

    cached_path = Path(gtf_path) if gtf_path is not None else None
    cache_root = user_data_dir(cache_dir) / "cache"
    clean_cache_files(
        cache_root,
        clean_cache,
        current_path=cached_path if cached_path and cached_path.parent == cache_root else None,
    )

    return output_bed


def _format_transcript_bed(
    record: _TranscriptRecord,
    selected_types: Optional[set],
) -> str:
    if selected_types is not None and record.gene_type not in selected_types:
        return ""
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
        record.gene_type,
    ]
    return "{}\n".format("\t".join(fields))


def regenerate_gencode_beds(
    output_dir: PathLike,
    specs: Iterable[Tuple[str, str, PathLike]],
    overwrite: bool = True,
    processes: int = 1,
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
                overwrite=overwrite,
                processes=processes,
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
