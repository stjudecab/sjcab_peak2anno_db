"""Install optional external BED resources into the user data cache."""

from __future__ import annotations

import gzip
import shutil
from pathlib import Path
from typing import Iterable, Optional, Tuple, Union

from ._download import ProgressCallback, download_file
from ._download_log import record_download_url
from ._registry import UnknownResourceError, data_root, user_data_dir

BLACKLIST_VERSION = "20230411"
BLACKLIST_DIR_NAME = "blacklists"
CGI_DIR_NAME = "cgi"
CGI_SPECIES = ("hg38", "hg19", "mm10", "mm9", "mm39")
CGI_TABLE = "cpgIslandExt"
CGI_URL_TEMPLATE = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/{species}/database/"
    + CGI_TABLE
    + ".txt.gz"
)


def blacklist_dir(data_dir: Optional[object] = None) -> Path:
    """Return the cache directory for blacklist BED files."""

    return user_data_dir(data_dir) / BLACKLIST_DIR_NAME


def cgi_dir(data_dir: Optional[object] = None) -> Path:
    """Return the cache directory for CpG island BED files."""

    return user_data_dir(data_dir) / CGI_DIR_NAME


def blacklist_path(name: str, data_dir: Optional[object] = None) -> Path:
    """Return a cache path for a blacklist by species or file name."""

    filename = name if name.endswith(".bed") else "{}-blacklist.bed".format(name)
    return blacklist_dir(data_dir) / filename


def cgi_path(species: str, data_dir: Optional[object] = None) -> Path:
    """Return the cache path for one species CpG island BED."""

    species = _normalize_cgi_species((species,))[0]
    return cgi_dir(data_dir) / "{}_cgi.bed".format(species)


def iter_blacklists() -> Tuple[str, ...]:
    """Return bundled blacklist file names."""

    return tuple(source.name for source in _blacklist_sources())


def install_blacklists(
    data_dir: Optional[object] = None,
    overwrite: bool = True,
) -> Path:
    """Install bundled blacklists into the cache with dated real files.

    Each bundled ``*.bed`` file is copied to ``*.bed.20230411``. The current
    ``*.bed`` name is refreshed as a relative symlink to the dated file, falling
    back to a copy on filesystems that do not support symlinks.
    """

    target_dir = blacklist_dir(data_dir)
    target_dir.mkdir(parents=True, exist_ok=True)

    for source in _blacklist_sources():
        versioned_path = target_dir / "{}.{}".format(source.name, BLACKLIST_VERSION)
        current_path = target_dir / source.name

        if overwrite or not versioned_path.exists():
            if versioned_path.exists() or versioned_path.is_symlink():
                versioned_path.unlink()
            shutil.copyfile(source, versioned_path)

        _refresh_link(current_path, versioned_path.name, versioned_path)

    return target_dir


def install_cgi(
    data_dir: Optional[object] = None,
    overwrite: bool = True,
) -> Path:
    """Install packaged CGI BED files into the cache."""

    target_dir = cgi_dir(data_dir)
    target_dir.mkdir(parents=True, exist_ok=True)

    for species in CGI_SPECIES:
        source = data_root() / CGI_DIR_NAME / "{}_cgi.bed".format(species)
        destination = target_dir / "{}_cgi.bed".format(species)
        if destination.exists() and not overwrite:
            continue
        shutil.copyfile(source, destination)

    return target_dir


def download_cgi(
    data_dir: Optional[object] = None,
    species: Optional[Union[str, Iterable[str]]] = None,
    overwrite: bool = True,
    progress: Optional[ProgressCallback] = None,
) -> Path:
    """Download UCSC CpG island tables and write BED-like files.

    UCSC ``cpgIslandExt`` rows are BED-like except for the leading ``bin``
    column. The generated ``*_cgi.bed`` files keep all columns after ``bin``.
    """

    selected_species = _normalize_cgi_species(species)
    target_dir = cgi_dir(data_dir)
    target_dir.mkdir(parents=True, exist_ok=True)

    for genome in selected_species:
        raw_path = target_dir / "UCSC_{}_{}.txt.gz".format(genome, CGI_TABLE)
        bed_path = target_dir / "{}_cgi.bed".format(genome)
        url = CGI_URL_TEMPLATE.format(species=genome)

        if overwrite or not raw_path.exists():
            if progress is None:
                _download_file(url, raw_path)
            else:
                _download_file(url, raw_path, progress=progress)
            record_download_url(url, data_dir=data_dir, destination=raw_path)

        if overwrite or not bed_path.exists():
            _write_cgi_bed(raw_path, bed_path)

    return target_dir


def _blacklist_sources() -> Tuple[Path, ...]:
    source_dir = data_root() / BLACKLIST_DIR_NAME
    if not source_dir.exists():
        return tuple()
    return tuple(sorted(source_dir.glob("*blacklist*.bed"), key=lambda path: path.name))


def _normalize_cgi_species(
    species: Optional[Union[str, Iterable[str]]]
) -> Tuple[str, ...]:
    if species is None:
        return CGI_SPECIES

    if isinstance(species, str):
        normalized = (species,)
    else:
        normalized = tuple(value for value in species)
    unknown = sorted(set(normalized) - set(CGI_SPECIES))
    if unknown:
        raise UnknownResourceError(
            "Unsupported CGI species {}. Supported species: {}".format(
                ", ".join(unknown), ", ".join(CGI_SPECIES)
            )
        )
    return normalized


def _download_file(
    url: str,
    destination: Path,
    progress: Optional[ProgressCallback] = None,
) -> None:
    download_file(url, destination, timeout=60, progress=progress)


def _write_cgi_bed(raw_path: Path, bed_path: Path) -> None:
    tmp_path = bed_path.with_name(bed_path.name + ".tmp")

    try:
        with gzip.open(raw_path, "rt", encoding="utf-8") as source, tmp_path.open(
            "w", encoding="utf-8"
        ) as output:
            for line_number, line in enumerate(source, start=1):
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue

                fields = line.split("\t")
                if len(fields) < 5:
                    raise ValueError(
                        "{}:{} is not a valid UCSC {} row".format(
                            raw_path, line_number, CGI_TABLE
                        )
                    )
                fields[4] = fields[4].replace(": ","_")
                output.write("{}\n".format("\t".join(fields[1:])))

        tmp_path.replace(bed_path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()


def _refresh_link(link_path: Path, target_name: str, fallback_source: Path) -> None:
    if link_path.exists() or link_path.is_symlink():
        link_path.unlink()

    try:
        link_path.symlink_to(target_name)
    except OSError:
        shutil.copyfile(fallback_source, link_path)
