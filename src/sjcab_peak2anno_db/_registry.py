"""Registry helpers for bundled and installed GeneBED resources."""

from __future__ import annotations

import os
import re
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Iterator, Mapping, Optional, TextIO, Tuple, Union

from ._config import configured_db_path

SUPPORTED_SPECIES = ("hg19", "hg38", "mm10", "mm9", "mm39", "sacCer3")
ANNOTATION_TYPES = ("gene",)
ISOFORM_SETS = ("all", "deduplong")
DATA_PATH_ENV_VAR = "SJCAB_PEAK2ANNO_DB_PATH"
DEFAULT_DATA_DIR = Path.home() / ".sjcab_peak2anno_db"

_SPECIES_ALIASES = {species.lower(): species for species in SUPPORTED_SPECIES}
_ANNOTATION_ALIASES = {
    "gene": "gene",
    "genes": "gene",
}
_ISOFORM_SET_ALIASES = {
    "all": "all",
    "full": "all",
    "all-isoforms": "all",
    "all_isoforms": "all",
    "deduplong": "deduplong",
    "dedup-long": "deduplong",
    "dedup_long": "deduplong",
    "longest": "deduplong",
    "longest-isoform": "deduplong",
    "longest_isoform": "deduplong",
}

_GENE_RESOURCES = (
    ("hg19", "v31lift37", "gencode.v31lift37.hg19.gene.bed.withtype"),
    ("hg38", "v31", "gencode.v31.hg38.gene.bed.withtype"),
    ("mm10", "vM22", "gencode.vM22.mm10.gene.bed.withtype"),
    ("mm9", "vM17", "gencode.vM17.mm9.gene.bed.withtype"),
    ("mm39", "vM39", "gencode.vM39.mm39.gene.bed.withtype"),
    ("sacCer3", "R64-1-1", "sacCer3.geneNames.bed"),
)

_DEFAULT_GENE_VERSIONS = {
    "hg19": "v31lift37",
    "hg38": "v31",
    "mm10": "vM22",
    "mm9": "vM17",
    "mm39": "vM39",
    "sacCer3": "R64-1-1",
}

PathLike = Union[str, os.PathLike]


class UnknownResourceError(KeyError):
    """Raised when a species, annotation type, or version is not available."""


@dataclass(frozen=True)
class AnnotationResource:
    """Description of one available annotation resource."""

    species: str
    isoform_set: str
    annotation: str
    version: str
    source_relative_path: str

    @property
    def is_derived(self) -> bool:
        """Whether the annotation is generated from a bundled GeneBED."""

        return self.isoform_set != "all" or self.annotation != "gene"

    @property
    def package_path(self) -> str:
        """Return the bundled gene path relative to the importable package."""

        return "data/{}".format(self.source_relative_path)

    @property
    def installed_relative_path(self) -> str:
        """Return the cache path relative to the selected user data directory."""

        return "genebed/{}/{}/{}.{}.bed".format(
            self.species, self.version, self.isoform_set, self.annotation
        )

    @property
    def relative_path(self) -> str:
        """Backward-compatible alias for the installed relative path."""

        return self.installed_relative_path

    @property
    def is_default(self) -> bool:
        """Whether this entry is the default for its species/type pair."""

        return self.version == default_version(
            self.species, self.annotation, self.isoform_set
        )


def data_root() -> Path:
    """Return the package directory containing bundled gene data files."""

    return Path(__file__).resolve().parent / "data"


def user_data_dir(data_dir: Optional[PathLike] = None) -> Path:
    """Return the user data directory used for generated annotations.

    Explicit ``data_dir`` wins. Otherwise ``SJCAB_PEAK2ANNO_DB_PATH`` is used
    when defined, falling back to ``~/.sjcab_peak2anno_db``.
    """

    if data_dir is not None:
        return Path(data_dir).expanduser()

    env_path = os.environ.get(DATA_PATH_ENV_VAR) or configured_db_path()
    if env_path:
        return Path(env_path).expanduser()
    return DEFAULT_DATA_DIR


def supported_species() -> Tuple[str, ...]:
    """Return supported species names."""

    return SUPPORTED_SPECIES


def iter_resources() -> Tuple[AnnotationResource, ...]:
    """Return all gene-backed installed resources."""

    resources_by_key = []
    for species, version, source_relative_path in _sorted_gene_resources():
        for isoform_set in ISOFORM_SETS:
            for annotation in ANNOTATION_TYPES:
                resources_by_key.append(
                    AnnotationResource(
                        species=species,
                        isoform_set=isoform_set,
                        annotation=annotation,
                        version=version,
                        source_relative_path=source_relative_path,
                    )
                )
    return tuple(resources_by_key)


def available() -> Mapping[str, Mapping[str, Mapping[str, Tuple[str, ...]]]]:
    """Return available versions grouped by species and annotation type."""

    return {
        species: {
            isoform_set: {
                annotation: versions(species, annotation, isoform_set)
                for annotation in ANNOTATION_TYPES
            }
            for isoform_set in ISOFORM_SETS
        }
        for species in SUPPORTED_SPECIES
    }


def versions(
    species: str, annotation: str = "gene", isoform_set: str = "all"
) -> Tuple[str, ...]:
    """Return available versions for a species and annotation type."""

    species, _annotation, _isoform_set = _normalize_resource_request(
        species, annotation, isoform_set
    )
    values = tuple(
        version
        for gene_species, version, _source in _sorted_gene_resources()
        if gene_species == species
    )
    if not values:
        raise UnknownResourceError("No versions available for {}.".format(species))
    return values


def default_version(
    species: str, annotation: str = "gene", isoform_set: str = "all"
) -> str:
    """Return the default version for a species and annotation type."""

    species, annotation, isoform_set = _normalize_resource_request(
        species, annotation, isoform_set
    )
    try:
        selected_version = _DEFAULT_GENE_VERSIONS[species]
    except KeyError as exc:
        raise UnknownResourceError(
            "No default version configured for {}.".format(species)
        ) from exc
    if selected_version not in versions(species, annotation, isoform_set):
        raise UnknownResourceError(
            "Default version {} is not available for {}/{}/{}.".format(
                selected_version, species, isoform_set, annotation
            )
        )
    return selected_version


def resource(
    species: str,
    annotation: str = "gene",
    version: Optional[str] = "default",
    isoform_set: str = "all",
) -> Path:
    """Return the bundled gene file path backing an annotation."""

    entry = _resolve(species, annotation, version, isoform_set)
    return data_root() / entry.source_relative_path


def installed_path(
    species: str,
    annotation: str = "gene",
    version: Optional[str] = "default",
    data_dir: Optional[PathLike] = None,
    isoform_set: str = "all",
) -> Path:
    """Return the expected path for an installed/generated annotation file."""

    entry = _resolve(species, annotation, version, isoform_set)
    selected_version = entry.version
    if version is None or version in {"default", "def", "latest"}:
        selected_version = "def"
    return (
        user_data_dir(data_dir)
        / "genebed"
        / entry.species
        / selected_version
        / "{}.{}.bed".format(entry.isoform_set, entry.annotation)
    )


def path(
    species: str,
    annotation: str = "gene",
    version: Optional[str] = "default",
    data_dir: Optional[PathLike] = None,
    prefer_user: bool = True,
    isoform_set: str = "all",
) -> Path:
    """Return a filesystem path for an annotation.

    This returns the installed user copy when present, otherwise the bundled
    package file. Use :func:`write_tss` and :func:`write_tes` to derive site
            annotations from a GeneBED when needed.
    """

    entry = _resolve(species, annotation, version, isoform_set)

    if prefer_user:
        cached = installed_path(
            entry.species,
            entry.annotation,
            version,
            data_dir,
            isoform_set=entry.isoform_set,
        )
        if cached.exists():
            return cached

    if entry.isoform_set == "all" and entry.annotation == "gene":
        return resource(entry.species, entry.annotation, entry.version)

    cached = installed_path(
        entry.species,
        entry.annotation,
        version,
        data_dir,
        isoform_set=entry.isoform_set,
    )
    raise UnknownResourceError(
        "{}/{}/{}/{} has not been generated at {}. Run "
        "sjcab_peak2anno_db.install_data() or `sjcab-peak2anno-db install` first.".format(
            entry.species,
            entry.isoform_set,
            entry.annotation,
            entry.version,
            cached,
        )
    )


@contextmanager
def as_file(
    species: str,
    annotation: str = "gene",
    version: Optional[str] = "default",
    isoform_set: str = "all",
) -> Iterator[Path]:
    """Yield the bundled gene file path backing an annotation."""

    entry = _resolve(species, annotation, version, isoform_set)
    yield resource(entry.species, "gene", entry.version)


def open_text(
    species: str,
    annotation: str = "gene",
    version: Optional[str] = "default",
    encoding: str = "utf-8",
    isoform_set: str = "all",
) -> TextIO:
    """Open the bundled gene file backing an annotation for text reading."""

    return resource(species, annotation, version, isoform_set).open(
        "r", encoding=encoding
    )


def open_binary(
    species: str,
    annotation: str = "gene",
    version: Optional[str] = "default",
    isoform_set: str = "all",
) -> BinaryIO:
    """Open the bundled gene file backing an annotation for binary reading."""

    return resource(species, annotation, version, isoform_set).open("rb")


def _resolve(
    species: str,
    annotation: str,
    version: Optional[str],
    isoform_set: str = "all",
) -> AnnotationResource:
    species, annotation, isoform_set = _normalize_resource_request(
        species, annotation, isoform_set
    )
    if version is None or version in {"default", "def", "latest"}:
        version = default_version(species, annotation, isoform_set)

    for gene_species, gene_version, source_relative_path in _sorted_gene_resources():
        if gene_species == species and gene_version == version:
            return AnnotationResource(
                species=gene_species,
                isoform_set=isoform_set,
                annotation=annotation,
                version=gene_version,
                source_relative_path=source_relative_path,
            )

    raise UnknownResourceError(
        "No resource available for {}/{}/{}/{}. Available versions: {}".format(
            species,
            isoform_set,
            annotation,
            version,
            ", ".join(versions(species, annotation, isoform_set)),
        )
    )


def _normalize_species(species: str) -> str:
    try:
        return _SPECIES_ALIASES[species.lower()]
    except KeyError as exc:
        raise UnknownResourceError(
            "Unsupported species {!r}. Supported species: {}".format(
                species, ", ".join(SUPPORTED_SPECIES)
            )
        ) from exc


def _normalize_annotation(annotation: str) -> str:
    key = annotation.lower()
    try:
        return _ANNOTATION_ALIASES[key]
    except KeyError as exc:
        raise UnknownResourceError(
            "Unsupported annotation type {!r}. Supported types: {}".format(
                annotation, ", ".join(ANNOTATION_TYPES)
            )
        ) from exc


def _normalize_isoform_set(isoform_set: str) -> str:
    key = isoform_set.lower()
    try:
        return _ISOFORM_SET_ALIASES[key]
    except KeyError as exc:
        raise UnknownResourceError(
            "Unsupported isoform set {!r}. Supported isoform sets: {}".format(
                isoform_set, ", ".join(ISOFORM_SETS)
            )
        ) from exc


def _normalize_resource_request(
    species: str, annotation: str, isoform_set: str = "all"
) -> Tuple[str, str, str]:
    species = _normalize_species(species)
    annotation_key = annotation.lower()
    isoform_set = _normalize_isoform_set(isoform_set)

    if annotation_key in _ISOFORM_SET_ALIASES:
        if isoform_set != "all" and _normalize_isoform_set(annotation_key) != isoform_set:
            raise UnknownResourceError(
                "Annotation {!r} conflicts with isoform set {!r}.".format(
                    annotation, isoform_set
                )
            )
        return species, "gene", _normalize_isoform_set(annotation_key)

    return species, _normalize_annotation(annotation), isoform_set


def _sorted_gene_resources() -> Tuple[Tuple[str, str, str], ...]:
    return tuple(
        sorted(
            _GENE_RESOURCES,
            key=lambda item: (
                SUPPORTED_SPECIES.index(item[0]),
                _version_key(item[1]),
                item[1],
            ),
        )
    )


def _version_key(version: str) -> Tuple[object, ...]:
    release = re.fullmatch(r"vM?(?P<major>\d+)(?:lift(?P<lift>\d+))?", version)
    if release:
        return (
            4,
            (int(release.group("major")), int(release.group("lift") or 0)),
            version.lower(),
        )

    yeast_release = re.match(
        r"R(?P<major>\d+)(?:-(?P<minor>\d+))?(?:-(?P<patch>\d+))?",
        version,
    )
    if yeast_release:
        return (
            4,
            (
                int(yeast_release.group("major")),
                int(yeast_release.group("minor") or 0),
                int(yeast_release.group("patch") or 0),
            ),
            version.lower(),
        )

    numbers = tuple(int(value) for value in re.findall(r"\d+", version))
    return (0, numbers, version.lower())
