"""Install bundled gene resources and derived annotations into a user cache."""

from __future__ import annotations

import json
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Optional, Tuple, Union

from ._download import ProgressCallback, report_progress
from ._config import load_config, parse_species_version_values
from ._derive import write_deduplong
from ._gencode import data_species_name, download_gencode_gtf_batch
from ._external import (
    BLACKLIST_VERSION,
    install_blacklists,
    install_cgi,
    iter_blacklists,
    CGI_SPECIES,
)
from ._registry import (
    ANNOTATION_TYPES,
    ISOFORM_SETS,
    AnnotationResource,
    DATA_PATH_ENV_VAR,
    default_version,
    iter_resources,
    resource,
    user_data_dir,
)
from ._regions import (
    GENCODE_FEATURE_LIST_ORDER,
    download_gencode_feature,
    gencode_feature_prefix,
)
from ._version import __version__

INSTALL_COMPONENTS = ("genebed", "feature", "blacklists", "cgi")
GENCODE_FEATURE_DIR_NAME = "feature"
DEFAULT_GENCODE_FEATURE_SPECS = (
    ("hg38", "v31"),
    ("hg19", "v31lift37"),
    ("mm10", "vM22"),
    ("mm39", "vM39"),
)


def install_data(
    data_dir: Optional[object] = None,
    overwrite: bool = True,
    components: Optional[Union[str, Iterable[str]]] = None,
    progress: Optional[ProgressCallback] = None,
    cache_dir: Optional[object] = None,
    ucsc_annotation: str = "ens",
    clean_cache: Optional[Union[bool, int]] = None,
    sizes_clean: Optional[bool] = None,
    processes: int = 4,
) -> Path:
    """Install selected resources into the user data directory.

    By default this installs the ``genebed`` (GeneBEDs) and ``feature``
    ``blacklists``, and ``cgi``. ``feature`` installs the bundled gene
    GeneBED files plus deduplicated longest-isoform GeneBEDs.
    ``cache_dir`` controls downloaded GTF placement; ``ucsc_annotation`` and
    ``clean_cache`` are forwarded to GENCODE feature generation.
    """

    if processes < 1:
        raise ValueError("processes must be at least 1")
    config = load_config()
    selected = _normalize_install_components(
        components
        if components is not None
        else (config.install_components if config.components_configured else None)
    )
    target_root = user_data_dir(data_dir)
    target_root.mkdir(parents=True, exist_ok=True)
    report_progress(progress, "install: started")

    if "genebed" in selected:
        report_progress(progress, "install: installing GeneBEDs")
        install_gencode_beds(target_root, overwrite=overwrite)
        report_progress(progress, "install: GeneBEDs done")
    if "feature" in selected:
        feature_kwargs = dict(
            overwrite=overwrite,
            progress=progress,
            cache_dir=target_root if cache_dir is None else cache_dir,
            ucsc_annotation=ucsc_annotation,
            clean_cache=clean_cache,
            sizes_clean=sizes_clean,
            processes=processes,
        )
        if config.feature_specs_configured:
            configured_specs = tuple(
                "{}:{}".format(feature_species, feature_version)
                for feature_species, feature_version in config.feature_specs
            )
            install_gencode_features(
                target_root,
                species=configured_specs,
                **feature_kwargs
            )
        else:
            install_gencode_features(target_root, **feature_kwargs)

    if "blacklists" in selected:
        report_progress(progress, "install: installing blacklists")
        install_blacklists(target_root, overwrite=overwrite)
        report_progress(progress, "install: blacklists done")
    if "cgi" in selected:
        report_progress(progress, "install: installing cgi")
        install_cgi(target_root, overwrite=overwrite)
        report_progress(progress, "install: cgi done")

    _write_manifest(target_root)
    report_progress(progress, "install: done")
    return target_root


def install_gencode_beds(
    data_dir: Optional[object] = None,
    overwrite: bool = True,
    species: Optional[str] = None,
    custom_name: Optional[str] = None,
) -> Path:
    """Install bundled GENCODE BEDs and deduplicated gene annotations."""

    target_root = user_data_dir(data_dir)
    target_root.mkdir(parents=True, exist_ok=True)
    if custom_name and not species:
        raise ValueError("A species/build is required when using custom_name.")

    gene_entries = [
        entry
        for entry in iter_resources()
        if entry.isoform_set == "all" and entry.annotation == "gene"
    ]
    if species and species.lower() not in {"all", "empty", "default"}:
        gene_entries = [
            entry for entry in gene_entries if entry.species.lower() == species.lower()
        ]
    if not gene_entries:
        raise ValueError("No bundled GeneBED resource found for {!r}.".format(species))
    for entry in gene_entries:
        storage_species = custom_name or entry.species
        _install_gene(
            entry,
            target_root,
            overwrite=overwrite,
            storage_species=storage_species,
        )
        if custom_name:
            _write_custom_name_mapping(target_root, entry.species, custom_name)

    for entry in gene_entries:
        _install_deduplong_gene(
            entry,
            target_root,
            overwrite=overwrite,
            storage_species=custom_name or entry.species,
        )

    selected_species = {entry.species for entry in gene_entries}
    for entry in iter_resources():
        if entry.species not in selected_species:
            continue
        if entry.version == default_version(
            entry.species, entry.annotation, entry.isoform_set
        ):
            _refresh_default(entry, target_root, custom_name or entry.species)

    return target_root


def install_gencode_features(
    data_dir: Optional[object] = None,
    overwrite: bool = True,
    species: Optional[Union[str, Iterable[str]]] = None,
    version: Optional[str] = None,
    output_dir: Optional[object] = None,
    prefix: Optional[str] = None,
    gtf_path: Optional[object] = None,
    gtf_url: Optional[str] = None,
    gene_bed: Optional[object] = None,
    promoter_bp: Union[int, str] = 2000,
    distal_bp: Union[int, str] = 50000,
    tes_bp: Union[int, str] = 2000,
    split_tss: bool = True,
    progress: Optional[ProgressCallback] = None,
    cache_dir: Optional[object] = None,
    ucsc_annotation: str = "ens",
    clean_cache: Optional[Union[bool, int]] = None,
    sizes_clean: Optional[bool] = None,
    custom_name: Optional[str] = None,
    collector_backend: str = "python",
    processes: int = 4,
    promoter_down_bp: Optional[Union[int, str]] = None,
    distal_down_bp: Optional[Union[int, str]] = None,
    tes_up_bp: Optional[Union[int, str]] = None,
) -> Path:
    """Install bundled derived annotations and downloaded GENCODE features.

    ``deduplong/gene`` keeps one longest isoform per gene. Downloaded
    feature-region BEDs from
    :func:`download_gencode_feature` are installed under
    ``feature/{species}/{version}/{prefix}`` with ``feature/{species}/def``
    pointing to the selected version/prefix directory.
    """

    if processes < 1:
        raise ValueError("processes must be at least 1")
    report_progress(progress, "install-feature: started")
    report_progress(progress, "install-feature: installing GENCODE BEDs")
    target_root = install_gencode_beds(data_dir, overwrite=overwrite)
    effective_cache_dir = target_root if cache_dir is None else cache_dir
    report_progress(progress, "install-feature: GENCODE BEDs done")

    specs = _feature_install_specs(species, version)
    if gtf_path is None:
        gtf_paths = download_gencode_gtf_batch(
            specs,
            output_dir=output_dir or target_root,
            gtf_url=gtf_url,
            log_data_dir=target_root,
            progress=progress,
            cache_dir=effective_cache_dir,
            ucsc_annotation=ucsc_annotation,
        )
    else:
        gtf_paths = tuple(Path(gtf_path).expanduser() for _ in specs)

    jobs = []
    for (feature_species, feature_version), selected_gtf_path in zip(specs, gtf_paths):
        report_progress(
            progress,
            "install-feature: installing feature set {} {}".format(
                feature_species, feature_version
            ),
        )
        jobs.append(
            {
                "species": feature_species,
                "version": feature_version,
                "data_dir": target_root,
                "overwrite": overwrite,
                "prefix": prefix,
                "source_dir": output_dir,
                "gtf_path": selected_gtf_path,
                "gtf_url": gtf_url,
                "gene_bed": gene_bed,
                "promoter_bp": promoter_bp,
                "distal_bp": distal_bp,
                "tes_bp": tes_bp,
                "split_tss": split_tss,
                "cache_dir": effective_cache_dir,
                "ucsc_annotation": ucsc_annotation,
                "clean_cache": clean_cache,
                "sizes_clean": sizes_clean,
                "custom_name": custom_name,
                "collector_backend": collector_backend,
                "promoter_down_bp": promoter_down_bp,
                "distal_down_bp": distal_down_bp,
                "tes_up_bp": tes_up_bp,
                "skip_existing": _can_skip_existing_feature_install(
                    source_dir=output_dir,
                    gtf_path=gtf_path,
                    gtf_url=gtf_url,
                    gene_bed=gene_bed,
                ),
            }
        )

    if processes > 1 and len(jobs) > 1:
        with ProcessPoolExecutor(max_workers=processes) as executor:
            tuple(executor.map(_install_gencode_feature_set_worker, jobs))
    else:
        for job in jobs:
            _install_gencode_feature_set_worker(job)

    for feature_species, feature_version in specs:
        report_progress(
            progress,
            "install-feature: feature set {} {} done".format(
                feature_species, feature_version
            ),
        )

    report_progress(progress, "install-feature: done")
    return target_root


def _install_gencode_feature_set_worker(job):
    arguments = dict(job)
    arguments["progress"] = _install_worker_progress
    arguments["processes"] = 1
    return install_gencode_feature_set(**arguments)


def _install_worker_progress(message: str) -> None:
    print(
        "[{}] {}".format(time.strftime("%H:%M:%S"), message),
        file=sys.stderr,
        flush=True,
    )


def install_gencode_feature_set(
    species: str,
    version: str,
    data_dir: Optional[object] = None,
    overwrite: bool = True,
    source_dir: Optional[object] = None,
    prefix: Optional[str] = None,
    gtf_path: Optional[object] = None,
    gtf_url: Optional[str] = None,
    gene_bed: Optional[object] = None,
    promoter_bp: Union[int, str] = 2000,
    distal_bp: Union[int, str] = 50000,
    tes_bp: Union[int, str] = 2000,
    split_tss: bool = True,
    progress: Optional[ProgressCallback] = None,
    skip_existing: bool = False,
    cache_dir: Optional[object] = None,
    ucsc_annotation: str = "ens",
    clean_cache: Optional[Union[bool, int]] = None,
    sizes_clean: Optional[bool] = None,
    custom_name: Optional[str] = None,
    collector_backend: str = "python",
    processes: int = 4,
    promoter_down_bp: Optional[Union[int, str]] = None,
    distal_down_bp: Optional[Union[int, str]] = None,
    tes_up_bp: Optional[Union[int, str]] = None,
) -> Path:
    """Install one downloaded GENCODE feature set into the cache."""

    if processes < 1:
        raise ValueError("processes must be at least 1")
    target_root = user_data_dir(data_dir)
    assembly_name = data_species_name(species, version, cache_dir=cache_dir)
    storage_species = custom_name or assembly_name
    if custom_name:
        _write_custom_name_mapping(target_root, assembly_name, custom_name)
    label = gencode_feature_prefix(promoter_bp=promoter_bp, prefix=prefix)
    species_dir = target_root / GENCODE_FEATURE_DIR_NAME / storage_species
    version_dir = species_dir / version
    feature_dir = version_dir / label
    feature_dir.mkdir(parents=True, exist_ok=True)

    if skip_existing and _is_preprocessed_feature_dir(feature_dir, label):
        _refresh_default_feature_dir(
            species_dir, version, label, feature_dir, overwrite=overwrite
        )
        report_progress(
            progress,
            "install-feature: using existing feature set {}".format(
                feature_dir
            ),
        )
        return feature_dir

    selected_gene_bed = (
        Path(gene_bed).expanduser()
        if gene_bed is not None
        else _existing_installed_gene_bed(target_root, storage_species, version)
    )
    gene_bed_output = (
        None
        if selected_gene_bed is not None
        else _installed_gene_bed_path(target_root, storage_species, version)
    )

    if source_dir is not None:
        source_feature_dir = _find_preprocessed_feature_dir(
            Path(source_dir).expanduser(), storage_species, version, label
        )
        if source_feature_dir is None:
            source_feature_dir = (
                Path(source_dir).expanduser() / storage_species / version / label
            )
            source_feature_dir.mkdir(parents=True, exist_ok=True)
            download_gencode_feature(
                species,
                version,
                source_feature_dir,
                gtf_path=gtf_path,
                gtf_url=gtf_url,
                gene_bed=selected_gene_bed,
                promoter_bp=promoter_bp,
                distal_bp=distal_bp,
                prefix=label,
                split_tss=split_tss,
                tes_bp=tes_bp,
                overwrite=overwrite,
                log_data_dir=target_root,
                progress=progress,
                cache_dir=cache_dir,
                ucsc_annotation=ucsc_annotation,
                clean_cache=clean_cache,
                sizes_clean=sizes_clean,
                data_species=storage_species,
                collector_backend=collector_backend,
                processes=processes,
                promoter_down_bp=promoter_down_bp,
                distal_down_bp=distal_down_bp,
                tes_up_bp=tes_up_bp,
                gene_bed_output=gene_bed_output,
            )
        _copy_preprocessed_feature_dir(
            source_feature_dir,
            feature_dir,
            label,
            overwrite=overwrite,
        )
    else:
        download_gencode_feature(
            species,
            version,
            feature_dir,
            gtf_path=gtf_path,
            gtf_url=gtf_url,
            gene_bed=selected_gene_bed,
            promoter_bp=promoter_bp,
            distal_bp=distal_bp,
            prefix=label,
            split_tss=split_tss,
            tes_bp=tes_bp,
            overwrite=overwrite,
            log_data_dir=target_root,
            progress=progress,
            cache_dir=cache_dir,
            ucsc_annotation=ucsc_annotation,
            clean_cache=clean_cache,
            sizes_clean=sizes_clean,
            data_species=storage_species,
            collector_backend=collector_backend,
            processes=processes,
            promoter_down_bp=promoter_down_bp,
            distal_down_bp=distal_down_bp,
            tes_up_bp=tes_up_bp,
            gene_bed_output=gene_bed_output,
        )
    if selected_gene_bed is None:
        selected_gene_bed = _installed_gene_bed_path(target_root, storage_species, version)
    deduplong_gene_bed = selected_gene_bed.parent / "deduplong.gene.bed"
    if selected_gene_bed.exists() and (overwrite or not deduplong_gene_bed.exists()):
        write_deduplong(selected_gene_bed, deduplong_gene_bed)
    _refresh_default_feature_dir(
        species_dir, version, label, feature_dir, overwrite=overwrite
    )
    report_progress(
        progress,
        "install-feature: refreshed default feature link {}".format(
            species_dir / "def"
        ),
    )
    return feature_dir


def _write_custom_name_mapping(
    target_root: Path, assembly_name: str, custom_name: str
) -> None:
    """Record the relationship between an assembly and a custom storage name."""

    mapping_path = target_root / "custom.name.tsv"
    existing = {}
    if mapping_path.exists():
        for line in mapping_path.read_text(encoding="utf-8").splitlines()[1:]:
            fields = line.split("\t")
            if len(fields) >= 2:
                existing[fields[0]] = fields[1]
    existing[assembly_name] = custom_name
    temporary = mapping_path.with_name(mapping_path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        handle.write("assembly\tcustom_name\n")
        for assembly, name in sorted(existing.items()):
            handle.write("{}\t{}\n".format(assembly, name))
    temporary.replace(mapping_path)


def update_data(data_dir: Optional[object] = None) -> Path:
    """Overwrite the user data directory with current package-derived files."""

    return install_data(data_dir=data_dir, overwrite=True)


def _normalize_install_components(
    components: Optional[Union[str, Iterable[str]]]
) -> Tuple[str, ...]:
    if components is None:
        return INSTALL_COMPONENTS
    if isinstance(components, str):
        values = (components,)
    else:
        values = tuple(components)

    unknown = sorted(set(values) - set(INSTALL_COMPONENTS))
    if unknown:
        raise ValueError(
            "Unsupported install component {}. Supported components: {}".format(
                ", ".join(unknown), ", ".join(INSTALL_COMPONENTS)
            )
        )
    return tuple(values)


def _feature_install_specs(
    species: Optional[Union[str, Iterable[str]]],
    version: Optional[str],
) -> Tuple[Tuple[str, str], ...]:
    if species is None or (isinstance(species, str) and species.lower() == "all"):
        if version is not None and str(version).lower() != "all":
            raise ValueError("version requires a species for install-feature.")
        return DEFAULT_GENCODE_FEATURE_SPECS

    species_value = ",".join(species) if not isinstance(species, str) else species
    specs = []
    for value, selected_version in parse_species_version_values(species_value, version):
        if selected_version is None or str(selected_version).lower() == "all":
            selected_version = default_version(value, "gene", "all")
        specs.append((value, selected_version))
    return tuple(specs)


def _can_skip_existing_feature_install(
    source_dir: Optional[object],
    gtf_path: Optional[object],
    gtf_url: Optional[str],
    gene_bed: Optional[object],
) -> bool:
    return (
        source_dir is None
        and gtf_path is None
        and gtf_url is None
        and gene_bed is None
    )


def _find_preprocessed_feature_dir(
    source_root: Path,
    species: str,
    version: str,
    label: str,
) -> Optional[Path]:
    candidates = (
        source_root / species / version / label,
        source_root / species / version,
        source_root / version / label,
        source_root / label,
        source_root,
    )
    for candidate in candidates:
        if _is_preprocessed_feature_dir(candidate, label):
            return candidate
    return None


def _is_preprocessed_feature_dir(path: Path, label: str) -> bool:
    if not path.is_dir():
        return False
    return all((path / name).exists() for name in _feature_file_names(label))


def _feature_file_names(label: str) -> Tuple[str, ...]:
    return tuple(
        "{}.{}.bed".format(label, region_type)
        for region_type in GENCODE_FEATURE_LIST_ORDER
    ) + ("order.lst",)


def _copy_preprocessed_feature_dir(
    source_dir: Path,
    destination_dir: Path,
    label: str,
    overwrite: bool,
) -> None:
    destination_dir.mkdir(parents=True, exist_ok=True)
    try:
        if source_dir.resolve() == destination_dir.resolve():
            return
    except OSError:
        pass
    for name in _feature_file_names(label):
        source = source_dir / name
        destination = destination_dir / name
        if destination.exists() and not overwrite:
            continue
        shutil.copy2(source, destination)
    for name in _optional_feature_file_names(label):
        source = source_dir / name
        if not source.exists():
            continue
        destination = destination_dir / name
        if destination.exists() and not overwrite:
            continue
        shutil.copy2(source, destination)


def _optional_feature_file_names(label: str) -> Tuple[str, ...]:
    return ("{}.promoter.bed".format(label),)


def _installed_gene_bed_path(target_root: Path, species: str, version: str) -> Path:
    return target_root / "bed" / species / version / "all.gene.bed"


def _existing_installed_gene_bed(
    target_root: Path,
    species: str,
    version: str,
) -> Optional[Path]:
    gene_bed = _installed_gene_bed_path(target_root, species, version)
    if gene_bed.exists():
        return gene_bed
    return None


def _install_gene(
    entry: AnnotationResource,
    target_root: Path,
    overwrite: bool,
    storage_species: Optional[str] = None,
) -> Path:
    destination = (
        user_data_dir(target_root)
        / "bed"
        / (storage_species or entry.species)
        / entry.version
        / "all.gene.bed"
    )
    if destination.exists() and not overwrite:
        return destination

    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(resource(entry.species, "gene", entry.version), destination)
    return destination


def _install_deduplong_gene(
    gene_entry: AnnotationResource,
    target_root: Path,
    overwrite: bool,
    storage_species: Optional[str] = None,
) -> Path:
    storage_species = storage_species or gene_entry.species
    destination = (
        user_data_dir(target_root)
        / "bed"
        / storage_species
        / gene_entry.version
        / "deduplong.gene.bed"
    )
    if destination.exists() and not overwrite:
        return destination

    gene_path = (
        user_data_dir(target_root)
        / "bed"
        / storage_species
        / gene_entry.version
        / "all.gene.bed"
    )
    write_deduplong(gene_path, destination)
    return destination


def _refresh_default(
    entry: AnnotationResource, target_root: Path, storage_species: Optional[str] = None
) -> None:
    storage_species = storage_species or entry.species
    version_dir = user_data_dir(target_root) / "bed" / storage_species / entry.version
    default_dir = user_data_dir(target_root) / "bed" / storage_species / "def"
    _refresh_dir_link(default_dir, entry.version, version_dir, overwrite=True)


def _refresh_default_feature_dir(
    species_dir: Path,
    version: str,
    label: str,
    feature_dir: Path,
    overwrite: bool,
) -> None:
    default_dir = species_dir / "def"
    _refresh_dir_link(
        default_dir,
        "{}/{}".format(version, label),
        feature_dir,
        overwrite=overwrite,
    )


def _refresh_dir_link(
    link_path: Path,
    target_name: str,
    fallback_source: Path,
    overwrite: bool,
) -> None:
    link_path.parent.mkdir(parents=True, exist_ok=True)
    if link_path.is_symlink():
        link_path.unlink()
    elif link_path.exists():
        if not overwrite:
            return
        if link_path.is_dir():
            shutil.rmtree(link_path)
        else:
            link_path.unlink()
    try:
        link_path.symlink_to(target_name, target_is_directory=True)
    except OSError:
        shutil.copytree(fallback_source, link_path, dirs_exist_ok=True)


def _write_manifest(target_root: Path) -> None:
    manifest = {
        "package": "sjcab_peak2anno_db",
        "package_version": __version__,
        "installed_at": datetime.now(timezone.utc).isoformat(),
        "data_dir": str(target_root),
        "env_var": DATA_PATH_ENV_VAR,
        "annotations": ANNOTATION_TYPES,
        "isoform_sets": ISOFORM_SETS,
        "blacklist_version": BLACKLIST_VERSION,
        "blacklists": [
            {
                "path": "blacklists/{}".format(name),
                "versioned_path": "blacklists/{}.{}".format(
                    name, BLACKLIST_VERSION
                ),
            }
            for name in iter_blacklists()
        ],
        "cgi_species": CGI_SPECIES,
        "cgi": [
            {
                "path": "cgi/{}_cgi.bed".format(species),
            }
            for species in CGI_SPECIES
        ],
        "resources": [
            {
                "species": entry.species,
                "isoform_set": entry.isoform_set,
                "annotation": entry.annotation,
                "version": entry.version,
                "default": entry.is_default,
                "path": entry.installed_relative_path,
                "source_gene": entry.source_relative_path,
            }
            for entry in iter_resources()
        ],
    }
    (target_root / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
