"""Command-line interface for sjcab_peak2anno_db."""

from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path
from typing import Optional, Sequence

from ._chromhmm import CHROMHMM_GENOMES, chromhmm_root, download_chromhmm
from ._config import parse_species_version_values
from ._dedup import dedup_gencode_bed, filter_gencode_bed
from ._external import CGI_SPECIES, download_cgi, install_blacklists, install_cgi
from ._gencode import (
    download_and_convert_gencode_gtf,
    ensembl_gtf_url_candidates,
    gtf_url_is_available,
    resolve_gencode_gtf_url,
)
from ._install import (
    DEFAULT_GENCODE_FEATURE_SPECS,
    INSTALL_COMPONENTS,
    install_data,
    install_gencode_beds,
    install_gencode_features,
    update_data,
)
from ._regions import download_gencode_feature, gencode_feature_prefix
from ._registry import (
    ANNOTATION_TYPES,
    ISOFORM_SETS,
    SUPPORTED_SPECIES,
    UnknownResourceError,
    iter_resources,
    path,
)
from ._segway import (
    download_segway,
    install_segway,
    segway_root,
    write_segway_liftover_script,
)
from ._liftover import write_liftover_script

_ANNOTATION_CHOICES = ANNOTATION_TYPES + ("deduplong",)
_GENCODE_FEATURE_SPECIES = ("hg38", "hg19", "mm10", "mm9", "mm39")
_GENCODE_SPECIES_HELP = (
    "Supported GENCODE builds: hg38, hg19."
)
_GENCODE_VERSION_HELP = (
    "Version examples: v31, v31lift37, vM23. Human releases: "
    "https://www.gencodegenes.org/human/releases.html Mouse releases: "
    "https://www.gencodegenes.org/mouse/releases.html"
)
_GENCODE_HELP_EPILOG = "{}\n{}".format(_GENCODE_SPECIES_HELP, _GENCODE_VERSION_HELP)
_PACKAGED_BLACKLIST_SPECIES = frozenset(
    {"ce10", "ce11", "dm3", "dm6", "hg18", "hg19", "hg38", "mm9", "mm10", "mm39", "saccer3"}
)
_RESOURCE_SPECIES_HELP = (
    "Available directly: hg38,hg19,mm39,mm10,mm9. Other UCSC builds will "
    "provide a liftOver bash script from hg38; check "
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/ for available."
)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(prog="sjcab-peak2anno-db")
    subparsers = parser.add_subparsers(dest="command", required=True)

    list_parser = subparsers.add_parser("list", help="List indexed resources.")
    list_parser.add_argument("-s", "--species", choices=SUPPORTED_SPECIES)
    list_parser.add_argument("-a", "--annotation", choices=ANNOTATION_TYPES)
    list_parser.add_argument("-i", "--isoform-set", choices=ISOFORM_SETS)

    path_parser = subparsers.add_parser("path", help="Print an annotation path.")
    path_parser.add_argument("species", choices=SUPPORTED_SPECIES)
    path_parser.add_argument("annotation", choices=_ANNOTATION_CHOICES)
    path_parser.add_argument("version", nargs="?", default="default")
    path_parser.add_argument("-v", "--ver", dest="version_option")
    path_parser.add_argument("-d", "--db-path", dest="data_dir", help="Generated annotation directory.")
    path_parser.add_argument(
        "-i",
        "--isoform-set",
        choices=ISOFORM_SETS,
        default="all",
        help="Use all isoforms or one longest isoform per gene.",
    )
    path_parser.add_argument(
        "-I",
        "--install",
        action="store_true",
        help="Generate all annotations before printing the path.",
    )

    install_parser = subparsers.add_parser(
        "install",
        help="Install/generate all annotations into the configured user data directory.",
    )
    install_parser.add_argument(
        "components",
        nargs="*",
        default=None,
        choices=INSTALL_COMPONENTS,
        metavar="COMPONENT",
        help=(
            "Components to install. Supported: {}. Defaults to all supported "
            "components.".format(
                ", ".join(INSTALL_COMPONENTS)
            )
        ),
    )
    install_parser.add_argument(
        "-c",
        "--component",
        action="append",
        choices=INSTALL_COMPONENTS,
        dest="component_options",
        help="Component to install. Can be passed more than once.",
    )
    install_parser.add_argument("-d", "--db-path", dest="data_dir", help="Generated annotation directory.")
    install_parser.add_argument(
        "--ucsc-source",
        choices=("ens", "refseq"),
        default="ens",
        help="UCSC gene annotation when a UCSC build is selected.",
    )
    install_parser.add_argument(
        "--clean-cache",
        nargs="?", const=-1, type=int, default=None, metavar="DAYS",
        help="Clean cache files older than DAYS (default 90); negative means immediately.",
    )
    install_parser.add_argument(
        "-j", "--processes", type=int, default=1,
        help="Number of worker processes for GTF/GeneBED and FeatureBED generation.",
    )
    install_parser.add_argument(
        "--sizes-clean", nargs="?", const=1, type=int, default=None, metavar="0|1",
        help="Create .sizes.clean files (default); use --sizes-clean 0 to disable.",
    )
    install_overwrite_group = install_parser.add_mutually_exclusive_group()
    install_overwrite_group.add_argument(
        "-n",
        "--no-overwrite",
        action="store_false",
        dest="overwrite",
        default=False,
        help="Do not rewrite existing generated files. This is the default.",
    )
    install_overwrite_group.add_argument(
        "--overwrite",
        action="store_true",
        dest="overwrite",
        help="Rewrite existing generated files.",
    )

    update_parser = subparsers.add_parser(
        "update", help="Regenerate all annotations into the configured user data directory."
    )
    update_parser.add_argument("species", nargs="?", help="Species/build to update.")
    update_parser.add_argument("-d", "--db-path", dest="data_dir", help="Generated annotation directory.")

    blacklist_parser = subparsers.add_parser(
        "install-blacklists",
        help="Install blacklist BEDs into the configured user data directory.",
    )
    blacklist_parser.add_argument("species_positional", nargs="?", metavar="SPECIES")
    blacklist_parser.add_argument("-d", "--db-path", dest="data_dir", help="Generated annotation directory.")
    blacklist_parser.add_argument(
        "-s", "--species", nargs="+",
        help=_RESOURCE_SPECIES_HELP,
    )
    blacklist_parser.add_argument(
        "--yes-liftover", action="store_true",
        help="Accept generating a CrossMap script from hg38 for unsupported builds.",
    )
    blacklist_parser.add_argument(
        "-n",
        "--no-overwrite",
        action="store_true",
        help="Do not rewrite existing blacklist files.",
    )

    install_cgi_parser = subparsers.add_parser(
        "install-cgi",
        help="Install packaged CGI BED files into the configured user data directory.",
    )
    install_cgi_parser.add_argument("species_positional", nargs="?", metavar="SPECIES")
    install_cgi_parser.add_argument("-d", "--db-path", dest="data_dir", help="Generated annotation directory.")
    install_cgi_parser.add_argument(
        "-s", "--species", nargs="+",
        help=_RESOURCE_SPECIES_HELP,
    )
    install_cgi_parser.add_argument(
        "--yes-liftover", action="store_true",
        help="Accept generating a CrossMap script from hg38 for unsupported builds.",
    )
    install_cgi_parser.add_argument(
        "-n",
        "--no-overwrite",
        action="store_true",
        help="Do not rewrite existing CGI files.",
    )

    feature_install_parser = subparsers.add_parser(
        "install-feature",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Install bundled GENCODE BEDs and downloaded feature annotations.",
        epilog=_GENCODE_HELP_EPILOG,
    )
    feature_install_parser.add_argument(
        "species",
        nargs="?",
        help="Genome build to install, empty or all for Gencode default feature builds.",
    )
    feature_install_parser.add_argument(
        "version",
        nargs="?",
        help="GENCODE/UCSC/Ensemble version to install, def/default/current/empty for default versions.",
    )
    feature_install_parser.add_argument(
        "--all",
        action="store_true",
        help="Install the default feature set for hg38, hg19, mm10, and mm39.",
    )
    feature_install_parser.add_argument(
        "-s",
        "--species",
        dest="species_option",
        nargs="+",
        help="One or more genome builds, comma-separated or from a .lst/.list file.",
    )
    feature_install_parser.add_argument(
        "-v",
        "--ver",
        "--version",
        dest="version_option",
        help="GENCODE version to install.",
    )
    feature_install_parser.add_argument(
        "-d", "--db-path", dest="data_dir", help="Generated annotation directory."
    )
    feature_install_parser.add_argument(
        "-o",
        "--output-dir",
        help="Preprocessed/download staging directory to install from.",
    )
    _add_feature_generation_arguments(
        feature_install_parser,
        include_output_dir=False,
        default_overwrite=False,
        include_name=True,
    )

    bed_install_parser = subparsers.add_parser(
        "install-genebed",
        help="Install bundled GENCODE BEDs and derived TSS/TES files.",
    )
    bed_install_parser.add_argument("species", nargs="?", help="Species/build to install.")
    bed_install_parser.add_argument("-d", "--db-path", dest="data_dir", help="Generated annotation directory.")
    bed_overwrite_group = bed_install_parser.add_mutually_exclusive_group()
    bed_overwrite_group.add_argument(
        "-n",
        "--no-overwrite",
        action="store_false",
        dest="overwrite",
        default=False,
        help="Do not rewrite existing generated files. This is the default.",
    )
    bed_overwrite_group.add_argument(
        "--overwrite",
        action="store_true",
        dest="overwrite",
        help="Rewrite existing generated files.",
    )

    cgi_parser = subparsers.add_parser(
        "download-cgi",
        help="Download UCSC cpgIslandExt tables and write CGI BED files.",
    )
    cgi_parser.add_argument("species_positional", nargs="?", metavar="SPECIES")
    cgi_parser.add_argument("-d", "--db-path", dest="data_dir", help="Generated annotation directory.")
    cgi_parser.add_argument(
        "-s",
        "--species",
        dest="species_option",
        choices=CGI_SPECIES,
        nargs="+",
        help=_RESOURCE_SPECIES_HELP,
    )
    cgi_parser.add_argument(
        "--yes-liftover", action="store_true",
        help="For another build, download hg38 and write a CrossMap script.",
    )
    cgi_parser.add_argument(
        "-n",
        "--no-overwrite",
        action="store_true",
        help="Do not rewrite existing CGI files.",
    )

    gencode_parser = subparsers.add_parser(
        "download-genebed",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Download or reuse a GENCODE GTF and write organized BED files.",
        epilog=_GENCODE_HELP_EPILOG,
    )
    gencode_parser.add_argument(
        "species",
        help="Genome build, for example hg38 or hg19.",
    )
    gencode_parser.add_argument(
        "version",
        nargs="?",
        default=None,
        help="GENCODE version, for example v31, v31lift37, or vM23.",
    )
    gencode_parser.add_argument(
        "-v", "--ver", dest="version_option",
        help="GENCODE/UCSC/Ensembl version; overrides the positional VERSION.",
    )
    gencode_parser.add_argument(
        "-o",
        "--output-dir",
        default=".",
        help="Directory for generated BED files. GTFs use the database cache.",
    )
    gencode_parser.add_argument("-d", "--db-path", dest="data_dir", help="Database/cache directory.")
    gencode_parser.add_argument("-g", "--gtf-path", help="Use an existing local GTF.")
    gencode_parser.add_argument("-u", "--url", help="Override the default GTF URL.")
    gencode_parser.add_argument(
        "--ucsc-source",
        choices=("ens", "refseq"),
        default="ens",
        help="UCSC gene annotation when a UCSC build is selected.",
    )
    gencode_parser.add_argument(
        "--clean-cache",
        nargs="?", const=-1, type=int, default=None, metavar="DAYS",
        help="Clean cache files older than DAYS (default 90); negative means immediately.",
    )
    gencode_parser.add_argument(
        "-b",
        "--output-bed",
        help="Explicit single generated GeneBED path.",
    )
    gencode_parser.add_argument(
        "-t",
        "--gene-type",
        action="append",
        help="Keep only this GENCODE gene_type. Can be passed more than once.",
    )
    gencode_parser.add_argument(
        "-n",
        "--no-overwrite",
        action="store_true",
        help="Do not rewrite existing GTF or BED files.",
    )
    gencode_parser.add_argument(
        "-j", "--processes", type=int, default=1,
        help="Number of worker processes for GTF-to-GeneBED conversion.",
    )

    _add_gencode_feature_parser(
        subparsers,
        "download-feature",
        "Download/convert GENCODE and write merged FeatureBED classes.",
    )

    _add_dedup_filter_parser(
        subparsers,
        "dedup-bed",
        "Select one isoform per gene and fall back to longest when unselected.",
    )
    _add_dedup_filter_parser(
        subparsers,
        "filter-bed",
        "Select one isoform per gene and omit genes without selector support.",
    )

    _add_chromhmm_parser(
        subparsers,
        "download-chromhmm",
        "Download Roadmap ChromHMM dense BED files.",
    )
    _add_chromhmm_parser(
        subparsers,
        "install-chromhmm",
        "Download Roadmap ChromHMM dense BED files into the user data directory.",
    )
    _add_segway_parser(
        subparsers,
        "download-segway",
        "Download Segway encyclopedia hg19 BED files.",
    )
    _add_segway_parser(
        subparsers,
        "install-segway",
        "Download Segway encyclopedia hg19 BED files into the user data directory.",
    )

    args = parser.parse_args(argv)

    try:
        if args.command == "list":
            _print_list(args.species, args.annotation, args.isoform_set)
            return 0
        if args.command == "path":
            if args.install:
                install_data(data_dir=args.data_dir, overwrite=False)
            print(
                path(
                    args.species,
                    args.annotation,
                    args.version_option or args.version,
                    data_dir=args.data_dir,
                    isoform_set=args.isoform_set,
                )
            )
            return 0
        if args.command == "install":
            target = install_data(
                data_dir=args.data_dir,
                overwrite=args.overwrite,
                components=_install_components_from_args(args),
                progress=_stderr_progress,
                cache_dir=args.data_dir,
                ucsc_annotation=args.ucsc_source,
                clean_cache=args.clean_cache,
                sizes_clean=args.sizes_clean,
                processes=args.processes,
            )
            print(target)
            return 0
        if args.command == "update":
            target = update_data(data_dir=args.data_dir)
            print(target)
            return 0
        if args.command == "install-blacklists":
            species = tuple(
                value for value in ((args.species_positional,) if args.species_positional else ())
            )
            species += tuple(args.species or ())
            packaged = tuple(value.lower() for value in species if value.lower() in _PACKAGED_BLACKLIST_SPECIES)
            unsupported = tuple(value for value in species if value.lower() not in _PACKAGED_BLACKLIST_SPECIES)
            install_species = packaged or (("hg38",) if unsupported else None)
            blacklist_kwargs = {
                "data_dir": args.data_dir,
                "overwrite": not args.no_overwrite,
            }
            if install_species is not None:
                blacklist_kwargs["species"] = install_species
            target = install_blacklists(**blacklist_kwargs)
            print(target)
            _write_external_liftover_scripts(
                target, "blacklists", unsupported, args.yes_liftover
            )
            return 0
        if args.command == "install-cgi":
            species = tuple(
                value for value in ((args.species_positional,) if args.species_positional else ())
            )
            species += tuple(args.species or ())
            packaged = tuple(value.lower() for value in species if value.lower() in CGI_SPECIES)
            unsupported = tuple(value.lower() for value in species if value.lower() not in CGI_SPECIES)
            install_species = packaged or (("hg38",) if unsupported else None)
            cgi_kwargs = {
                "data_dir": args.data_dir,
                "overwrite": not args.no_overwrite,
            }
            if install_species is not None:
                cgi_kwargs["species"] = install_species
            target = install_cgi(**cgi_kwargs)
            print(target)
            _write_external_liftover_scripts(
                target, "cgi", unsupported, args.yes_liftover
            )
            return 0
        if args.command == "install-feature":
            feature_species, feature_version = _gencode_feature_install_scope(args)
            if args.dry_run:
                dry_specs = (
                    DEFAULT_GENCODE_FEATURE_SPECS
                    if feature_species is None
                    else _gencode_feature_specs(feature_species, feature_version)
                )
                for dry_species, dry_version in dry_specs:
                    print(_dry_run_gtf_url(dry_species, dry_version, args.data_dir))
                return 0
            target = install_gencode_features(
                data_dir=args.data_dir,
                overwrite=args.overwrite,
                species=feature_species,
                version=feature_version,
                output_dir=args.output_dir,
                gtf_path=args.gtf_path,
                gtf_url=args.url,
                gene_bed=args.gene_bed,
                promoter_bp=args.promoter_bp,
                distal_bp=args.distal_bp,
                prefix=args.prefix,
                split_tss=not args.include_tss_base,
                tes_bp=args.tes_bp,
                progress=_stderr_progress,
                custom_name=args.name,
                clean_cache=args.clean_cache,
                sizes_clean=args.sizes_clean,
                processes=args.processes,
            )
            print(target)
            return 0
        if args.command == "install-genebed":
            target = install_gencode_beds(
                data_dir=args.data_dir,
                overwrite=args.overwrite,
            )
            print(target)
            return 0
        if args.command == "download-cgi":
            species = tuple(
                value for value in ((args.species_positional,) if args.species_positional else ())
            )
            species += tuple(args.species_option or ())
            packaged = tuple(value.lower() for value in species if value.lower() in CGI_SPECIES)
            unsupported = tuple(value.lower() for value in species if value.lower() not in CGI_SPECIES)
            if unsupported and not args.yes_liftover:
                raise UnknownResourceError(
                    "CGI has no packaged files for {}. Re-run with --yes-liftover "
                    "to download hg38 and write CrossMap scripts.".format(
                        ", ".join(unsupported)
                    )
                )
            cgi_kwargs = {
                "data_dir": args.data_dir,
                "overwrite": not args.no_overwrite,
                "progress": _stderr_progress,
            }
            if species:
                cgi_kwargs["species"] = tuple(
                    dict.fromkeys(packaged + (("hg38",) if unsupported else ()))
                )
            target = download_cgi(**cgi_kwargs)
            print(target)
            _write_external_liftover_scripts(
                target, "cgi", unsupported, args.yes_liftover
            )
            return 0
        if args.command == "download-genebed":
            version = args.version_option or args.version or "def"
            target = download_and_convert_gencode_gtf(
                args.species,
                version,
                args.output_dir,
                gtf_path=args.gtf_path,
                gtf_url=args.url,
                output_bed=args.output_bed,
                gene_types=args.gene_type,
                overwrite=not args.no_overwrite,
                progress=_stderr_progress,
                clean_cache=args.clean_cache,
                processes=args.processes,
            )
            print(target)
            return 0
        if args.command == "download-feature":
            feature_species = args.species or args.species_option
            if not feature_species:
                raise ValueError("download-feature requires a species.")
            if isinstance(feature_species, (list, tuple)):
                feature_species = ",".join(feature_species)
            specs = _gencode_feature_specs(
                feature_species, args.version_option or args.version or "def"
            )
            if args.dry_run:
                for dry_species, dry_version in specs:
                    print(_dry_run_gtf_url(dry_species, dry_version, args.data_dir))
                return 0
            multi_spec = len(specs) > 1
            for feature_species, feature_version in specs:
                output_dir = _gencode_feature_output_dir(
                    args.output_dir,
                    feature_species,
                    feature_version,
                    args.promoter_bp,
                    args.prefix,
                    multi_spec,
                )
                targets = download_gencode_feature(
                    feature_species,
                    feature_version,
                    output_dir,
                    gtf_path=args.gtf_path,
                    gtf_url=args.url,
                    gene_bed=args.gene_bed,
                    promoter_bp=args.promoter_bp,
                    distal_bp=args.distal_bp,
                    prefix=args.prefix,
                    split_tss=not args.include_tss_base,
                    tes_bp=args.tes_bp,
                    overwrite=not args.no_overwrite,
                    progress=_stderr_progress,
                    cache_dir=args.data_dir,
                    ucsc_annotation=args.ucsc_source,
                    clean_cache=args.clean_cache,
                    sizes_clean=args.sizes_clean,
                    processes=args.processes,
                )
                for name in sorted(targets):
                    if multi_spec:
                        print(
                            "{}\t{}\t{}\t{}".format(
                                feature_species,
                                feature_version,
                                name,
                                targets[name],
                            )
                        )
                    else:
                        print("{}\t{}".format(name, targets[name]))
            return 0
        if args.command in {"dedup-bed", "filter-bed"}:
            function = (
                dedup_gencode_bed
                if args.command == "dedup-bed"
                else filter_gencode_bed
            )
            targets = function(
                args.method,
                args.selector,
                output_dir=args.output_dir,
                gene_bed=args.gene_bed,
                species=args.species,
                version=args.version_option or args.version,
                data_dir=args.data_dir,
                promoter_bp=args.promoter_bp,
                inclusive=not args.exclusive,
                output_prefix=args.prefix,
                gene_key=args.gene_key,
            )
            for name in sorted(targets):
                print("{}\t{}".format(name, targets[name]))
            return 0
        if args.command in {"download-chromhmm", "install-chromhmm"}:
            chromhmm_data_dir = (
                args.output_dir if args.command == "download-chromhmm" else args.data_dir
            )
            requested_genome = str(
                args.species_positional or args.genome
            ).lower()
            source_genome = requested_genome if requested_genome in CHROMHMM_GENOMES else "hg38"
            if source_genome != requested_genome and not args.yes_liftover:
                raise UnknownResourceError(
                    "ChromHMM has hg19/hg38 files. Re-run with --yes-liftover "
                    "to generate a CrossMap script for {}.".format(requested_genome)
                )
            targets = download_chromhmm(
                data_dir=chromhmm_data_dir,
                model=args.model,
                genome=source_genome,
                ids=args.ids,
                tissue=args.tissue,
                cellline=args.cellline,
                overwrite=not args.no_overwrite,
                progress=_stderr_progress,
            )
            for eid in sorted(targets):
                print("{}\t{}".format(eid, targets[eid]))
            if source_genome != requested_genome:
                root = chromhmm_root(chromhmm_data_dir)
                source_dir = root / source_genome / "{}state".format(args.model)
                output_dir = root / ".liftover" / requested_genome / "{}state".format(args.model)
                script = write_liftover_script(
                    source_dir, output_dir, source_genome, requested_genome,
                    root / "liftover_{}_to_{}.sh".format(source_genome, requested_genome),
                    install_root=(
                        root / requested_genome / "{}state".format(args.model)
                        if args.command == "install-chromhmm"
                        else None
                    ),
                    cache_dir=root.parent,
                )
                print("liftover_script\t{}".format(script))
            return 0
        if args.command in {"download-segway", "install-segway"}:
            segway_data_dir = (
                args.output_dir if args.command == "download-segway" else args.data_dir
            )
            requested_genome = args.species_positional or args.genome
            download_genome = requested_genome
            liftover_script = None
            if requested_genome.lower() != "hg19":
                if not args.yes_liftover and not _confirm_segway_liftover(
                    requested_genome
                ):
                    raise UnknownResourceError(
                        "Segway encyclopedia downloads are hg19 only. Re-run with "
                        "--yes-liftover to download hg19 and write a liftover "
                        "script for {}.".format(requested_genome)
                    )
                download_genome = "hg19"

            segway_function = (
                install_segway
                if args.command == "install-segway"
                else download_segway
            )
            targets = segway_function(
                data_dir=segway_data_dir,
                genome=download_genome,
                names=args.names,
                tissue=args.tissue,
                cellline=args.cellline,
                all_celltypes=args.all_celltypes,
                include_encyclopedia=args.include_encyclopedia,
                include_caas=args.include_caas,
                include_label_info=args.include_label_info,
                overwrite=args.overwrite,
                progress=_stderr_progress,
            )
            if requested_genome.lower() != "hg19":
                liftover_script = write_segway_liftover_script(
                    data_dir=segway_data_dir,
                    target_genome=requested_genome,
                    install_root=(
                        segway_root(segway_data_dir)
                        if args.command == "install-segway"
                        else None
                    ),
                )
            for name in sorted(targets):
                print("{}\t{}".format(name, targets[name]))
            if liftover_script is not None:
                print("liftover_script\t{}".format(liftover_script))
            return 0
    except (UnknownResourceError, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        return 2

    parser.error("Unsupported command: {}".format(args.command))
    return 2


def _add_gencode_feature_parser(
    subparsers: argparse._SubParsersAction,
    name: str,
    help_text: str,
    aliases: Sequence[str] = (),
) -> None:
    parser = subparsers.add_parser(
        name,
        aliases=aliases,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help=help_text,
        epilog=_GENCODE_HELP_EPILOG,
    )
    parser.add_argument(
        "species",
        nargs="?",
        help="Genome build, comma-separated or from a .lst/.list file.",
    )
    parser.add_argument(
        "version",
        nargs="?",
        default=None,
        help="GENCODE/UCSC/Ensemble version, for example v31, v31lift37, or vM23.",
    )
    parser.add_argument(
        "-v", "--ver", dest="version_option",
        help="GENCODE/UCSC/Ensembl version; overrides the positional VERSION.",
    )
    parser.add_argument(
        "-s", "--species", dest="species_option", nargs="+",
        help="One or more genome builds, comma-separated or from a .lst/.list file.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default=".",
        help="Directory for generated GeneBED and FeatureBED files.",
    )
    parser.add_argument("-d", "--db-path", dest="data_dir", help="Database/cache directory.")
    _add_feature_generation_arguments(parser, include_output_dir=False)


def _add_dedup_filter_parser(
    subparsers: argparse._SubParsersAction,
    name: str,
    help_text: str,
) -> None:
    parser = subparsers.add_parser(name, help=help_text)
    parser.add_argument("species_positional", nargs="?", metavar="SPECIES")
    parser.add_argument(
        "species",
        nargs="?",
        help="Genome build to resolve under the data directory when --gene-bed is absent.",
    )
    parser.add_argument(
        "version",
        nargs="?",
        default="default",
        help="GENCODE/UCSC/Ensemble version to resolve. Defaults to the species default.",
    )
    parser.add_argument(
        "-v", "--ver", dest="version_option",
        help="GENCODE/UCSC/Ensembl version; overrides the positional VERSION.",
    )
    parser.add_argument(
        "-b",
        "--gene-bed",
        help="Input all-isoform GeneBED. Overrides species/version lookup.",
    )
    method_help = "Selection method: longcol5, long, peak, isoID, isoexp, or perover."
    if name == "dedup-bed":
        parser.add_argument("-m", "--method", default="longcol5", help=method_help)
    else:
        parser.add_argument("-m", "--method", required=True, help=method_help)
    parser.add_argument(
        "-i",
        "--input",
        dest="selector",
        help="Selector file; not needed for longcol5 or long.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default=".",
        help="Directory for generated gene/tss/tes BED files.",
    )
    parser.add_argument(
        "-d",
        "--db-path",
        dest="data_dir",
        help="Generated annotation directory for species/version lookup.",
    )
    parser.add_argument(
        "-p",
        "--promoter-bp",
        default="2kb",
        help="Promoter half-window around TSS for peak/perover. Defaults to 2kb.",
    )
    parser.add_argument(
        "-P",
        "--prefix",
        help=(
            "Output prefix. Defaults to {species}.{version}.{dedup/filter}{method} "
            "or the input BED basename plus .{dedup/filter}{method}."
        ),
    )
    parser.add_argument(
        "-K",
        "--gene-key",
        choices=("symbol", "ensid"),
        default="symbol",
        help="Group isoforms by gene symbol or Ensembl/GENCODE gene ID. Defaults to symbol.",
    )
    parser.add_argument(
        "--exclusive",
        action="store_true",
        help=(
            "Use exact transcript ID matching for isoID/isoexp and require BED "
            "features to be fully contained in the promoter for peak/perover."
        ),
    )


def _add_feature_generation_arguments(
    parser: argparse.ArgumentParser,
    include_output_dir: bool,
    default_overwrite: bool = True,
    include_name: bool = False,
) -> None:
    if include_output_dir:
        parser.add_argument(
            "-o",
            "--output-dir",
            default=".",
        help="Directory for generated GeneBED and FeatureBED files.",
        )
    parser.add_argument("-g", "--gtf-path", help="Use an existing local GTF.")
    parser.add_argument("-u", "--url", help="Override the default GTF URL.")
    parser.add_argument(
        "--ucsc-source",
        choices=("ens", "refseq"),
        default="ens",
        help="UCSC gene annotation when a UCSC build is selected.",
    )
    parser.add_argument(
        "--clean-cache",
        nargs="?", const=-1, type=int, default=None, metavar="DAYS",
        help="Clean cache files older than DAYS (default 90); negative means immediately.",
    )
    parser.add_argument(
        "--sizes-clean", nargs="?", const=1, type=int, default=None, metavar="0|1",
        help="Create .sizes.clean files (default); use --sizes-clean 0 to disable.",
    )
    parser.add_argument(
        "-j", "--processes", type=int, default=1,
        help="Number of worker processes for GTF/GeneBED and FeatureBED generation.",
    )
    parser.add_argument(
        "-b",
        "--gene-bed",
        help=(
            "Use an existing GeneBED for TSS-flank files instead of converting "
            "one from the GTF."
        ),
    )
    parser.add_argument(
        "-p",
        "--promoter-bp",
        default="2kb",
        help="Promoter flank size, for example 2kb or 2000.",
    )
    parser.add_argument(
        "-D",
        "--distal-bp",
        default="50kb",
        help="Distal flank size, for example 50kb or 50000.",
    )
    parser.add_argument(
        "-e",
        "--tes-bp",
        default="2kb",
        help="TES flank size for tes.bed, for example 2kb or 2000.",
    )
    parser.add_argument(
        "-P",
        "--prefix",
        help="Output file prefix. Defaults to the promoter size label.",
    )
    parser.add_argument(
        "-B",
        "--include-tss-base",
        action="store_true",
        help="Start right-side regions at the TSS instead of after the TSS base.",
    )
    if include_name:
        parser.add_argument(
            "-name",
            dest="name",
            help="Customize the installed species name and record it in custom.name.tsv.",
        )
    parser.add_argument(
        "-dry-run",
        dest="dry_run",
        action="store_true",
        help="Resolve and print the downloadable GTF URL without downloading.",
    )
    if default_overwrite:
        parser.add_argument(
            "-n",
            "--no-overwrite",
            action="store_true",
            help="Do not rewrite existing GTF or BED files.",
        )
    else:
        overwrite_group = parser.add_mutually_exclusive_group()
        overwrite_group.add_argument(
            "-n",
            "--no-overwrite",
            action="store_false",
            dest="overwrite",
            default=False,
            help="Do not rewrite existing GTF or BED files. This is the default.",
        )
        overwrite_group.add_argument(
            "--overwrite",
            action="store_true",
            dest="overwrite",
            help="Rewrite existing GTF or BED files.",
        )


def _add_chromhmm_parser(
    subparsers: argparse._SubParsersAction,
    name: str,
    help_text: str,
) -> None:
    parser = subparsers.add_parser(
        name,
        help=help_text,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=_RESOURCE_SPECIES_HELP,
    )
    parser.add_argument("species_positional", nargs="?", metavar="SPECIES")
    if name == "download-chromhmm":
        parser.add_argument(
            "-o",
            "--output-dir",
            default=".",
            help="Directory for downloaded ChromHMM files.",
        )
    else:
        parser.add_argument(
            "-d", "--db-path", dest="data_dir", help="Generated annotation directory."
        )
    parser.add_argument(
        "-m",
        "--model",
        choices=("15", "18", "25"),
        default="18",
        help="Roadmap ChromHMM state model. Defaults to 18.",
    )
    parser.add_argument(
        "-s", "--species", dest="genome", metavar="SPECIES", default="hg19",
        help=_RESOURCE_SPECIES_HELP,
    )
    parser.add_argument(
        "--yes-liftover", action="store_true",
        help="Accept generating a CrossMap script from hg38 for another build.",
    )
    parser.add_argument(
        "-i",
        "--id",
        "--ids",
        dest="ids",
        action="append",
        help="Epigenome ID list, for example E001,E063. Can be repeated.",
    )
    parser.add_argument(
        "-t",
        "--tissue",
        action="append",
        help="Fuzzy tissue/group keyword, for example brain or blood.",
    )
    parser.add_argument(
        "-c",
        "--cellline",
        action="append",
        help="Fuzzy cell-line/sample keyword, for example H1 or GM12878.",
    )
    parser.add_argument(
        "-n",
        "--no-overwrite",
        action="store_true",
        help="Do not rewrite existing ChromHMM files.",
    )


def _add_segway_parser(
    subparsers: argparse._SubParsersAction,
    name: str,
    help_text: str,
) -> None:
    parser = subparsers.add_parser(
        name,
        help=help_text,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            _RESOURCE_SPECIES_HELP
        ),
    )
    parser.add_argument("species_positional", nargs="?", metavar="SPECIES")
    if name == "download-segway":
        parser.add_argument(
            "-o",
            "--output-dir",
            default=".",
            help="Directory for downloaded Segway files.",
        )
    else:
        parser.add_argument(
            "-d", "--db-path", dest="data_dir", help="Generated annotation directory."
        )
    parser.add_argument(
        "-s", "--species", dest="genome", metavar="SPECIES", default="hg19",
        help=_RESOURCE_SPECIES_HELP,
    )
    parser.add_argument(
        "-i",
        "--id",
        "--ids",
        "--sample",
        "--samples",
        dest="names",
        action="append",
        help=(
            "Segway sample/resource list, for example GM12878,H1-HESC,caas. "
            "Can be repeated."
        ),
    )
    parser.add_argument(
        "-t",
        "--tissue",
        action="append",
        help="Fuzzy tissue/sample keyword, for example brain or liver.",
    )
    parser.add_argument(
        "-c",
        "--cellline",
        action="append",
        help="Fuzzy cell-line/sample keyword, for example H1 or GM12878.",
    )
    parser.add_argument(
        "--all-celltypes",
        action="store_true",
        help="Download all discovered cell type-specific Segway BED files.",
    )
    parser.add_argument(
        "--include-encyclopedia",
        action="store_true",
        help="Download segway_encyclopedia.bed.gz.",
    )
    parser.add_argument(
        "--include-caas",
        action="store_true",
        help="Download caas.bed.gz.",
    )
    parser.add_argument(
        "--include-label-info",
        action="store_true",
        help="Download label_info.tab.",
    )
    parser.add_argument(
        "--yes-liftover",
        action="store_true",
        help=(
            "For non-hg19 --species requests, download hg19 files and write a "
            "CrossMap liftover script without prompting."
        ),
    )
    overwrite_group = parser.add_mutually_exclusive_group()
    overwrite_group.add_argument(
        "-n",
        "--no-overwrite",
        action="store_false",
        dest="overwrite",
        default=False,
        help="Do not rewrite existing Segway files. This is the default.",
    )
    overwrite_group.add_argument(
        "--overwrite",
        action="store_true",
        dest="overwrite",
        help="Rewrite existing Segway files.",
    )


def _confirm_segway_liftover(genome: str) -> bool:
    if not sys.stdin.isatty():
        return False
    prompt = (
        "Segway encyclopedia downloads are hg19 only. Download hg19 and write "
        "a CrossMap liftover script for {}? [y/N] ".format(genome)
    )
    response = input(prompt).strip().lower()
    return response in {"y", "yes"}


def _write_external_liftover_scripts(
    source_dir: Path,
    resource_name: str,
    target_genomes: Sequence[str],
    yes_liftover: bool,
) -> None:
    """Write hg38-to-target helpers for resources absent from the package."""

    for target_genome in target_genomes:
        if not yes_liftover and not _confirm_generic_liftover(
            resource_name, target_genome
        ):
            raise UnknownResourceError(
                "{} has no packaged {} file. Re-run with --yes-liftover to "
                "write a CrossMap script from hg38.".format(
                    target_genome, resource_name
                )
            )
        output_dir = source_dir / ".liftover" / target_genome
        script = write_liftover_script(
            source_dir,
            output_dir,
            "hg38",
            target_genome,
            source_dir / "liftover_hg38_to_{}.sh".format(target_genome),
            install_root=source_dir,
            rename_genome_prefix=resource_name in {"blacklists", "cgi"},
            cache_dir=source_dir.parent,
        )
        print("liftover_script\t{}".format(script))


def _confirm_generic_liftover(resource_name: str, genome: str) -> bool:
    if not sys.stdin.isatty():
        return False
    response = input(
        "{} has no packaged {} file. Generate a CrossMap script from hg38 "
        "for {}? [y/N] ".format(genome, resource_name, genome)
    )
    return response.strip().lower() in {"y", "yes"}


def _install_components_from_args(args: argparse.Namespace) -> Optional[Sequence[str]]:
    components = list(args.components or ())
    components.extend(args.component_options or ())
    return components or None


def _gencode_feature_install_scope(args: argparse.Namespace) -> tuple:
    if args.all:
        if args.species or args.version or args.species_option or args.version_option:
            raise ValueError(
                "install-feature --all cannot be combined with "
                "--species or --version."
            )
        return None, None
    species = args.species or args.species_option
    version = args.version or args.version_option
    if not species:
        raise ValueError(
            "install-feature requires species, or --all."
        )
    if version is None:
        species_text = ",".join(species) if isinstance(species, (list, tuple)) else species
        parsed = parse_species_version_values(species_text)
        version = None if any(item_version is not None for _, item_version in parsed) else "def"
    return species, version


def _dry_run_gtf_url(species: str, version: str, cache_dir: Optional[str]) -> str:
    """Validate the initial URL and interactively resolve 404 alternatives."""

    url = resolve_gencode_gtf_url(species, version, cache_dir=cache_dir)
    if gtf_url_is_available(url):
        return url

    candidates = []
    for fields, candidate_url in ensembl_gtf_url_candidates(
        species, version, cache_dir=cache_dir
    ):
        if gtf_url_is_available(candidate_url):
            candidates.append((fields, candidate_url))
    if not candidates:
        raise ValueError(
            "No non-404 GTF URL found after resolving Ensembl species catalog for {!r}.".format(
                species
            )
        )
    if len(candidates) == 1:
        return candidates[0][1]

    print("Multiple non-404 GTF URLs found; select one:", file=sys.stderr)
    print("#\tassembly\tspecies\tdivision\tname\tassembly_accession", file=sys.stderr)
    for index, (fields, candidate_url) in enumerate(candidates, start=1):
        print(
            "{}\t{}\t{}\t{}\t{}\t{}\n  {}".format(
                index, fields[3], fields[0], fields[2], fields[1], fields[4], candidate_url
            ),
            file=sys.stderr,
        )
    while True:
        choice = input("Select GTF URL [1-{}]: ".format(len(candidates))).strip()
        try:
            selected = int(choice)
        except ValueError:
            selected = 0
        if 1 <= selected <= len(candidates):
            return candidates[selected - 1][1]


def _gencode_feature_specs(species: str, version: Optional[str]) -> tuple:
    if species.lower() == "all" and (version or "all").lower() == "all":
        return DEFAULT_GENCODE_FEATURE_SPECS
    if species.lower() == "all":
        raise ValueError("download-feature all requires version all.")
    specs = parse_species_version_values(species, version)
    selected = []
    for feature_species, feature_version in specs:
        if feature_version is None or feature_version.lower() == "all":
            feature_version = _gencode_default_feature_version(feature_species)
        selected.append((feature_species, feature_version))
    return tuple(selected)


def _gencode_default_feature_version(species: str) -> str:
    for feature_species, feature_version in DEFAULT_GENCODE_FEATURE_SPECS:
        if feature_species == species:
            return feature_version
    raise ValueError(
        "No default GENCODE feature version configured for {}.".format(species)
    )


def _gencode_feature_output_dir(
    output_dir: str,
    species: str,
    version: str,
    promoter_bp: str,
    prefix: Optional[str],
    multi_spec: bool,
) -> str:
    if not multi_spec:
        return output_dir
    label = _gencode_feature_output_label(promoter_bp, prefix)
    return str(Path(output_dir).expanduser() / species / version / label)


def _gencode_feature_output_label(promoter_bp: str, prefix: Optional[str]) -> str:
    return gencode_feature_prefix(promoter_bp=promoter_bp, prefix=prefix)


def _stderr_progress(message: str) -> None:
    if message.startswith("download ") and message.endswith(":"):
        print(_log_time_prefix() + message, end="", file=sys.stderr, flush=True)
        return
    if re.fullmatch(r"\d+\.\.", message):
        print(message, end="", file=sys.stderr, flush=True)
        return
    if message == "done":
        print(message, file=sys.stderr, flush=True)
        return
    print(_log_time_prefix() + message, file=sys.stderr, flush=True)


def _log_time_prefix() -> str:
    """Return a second-resolution timestamp for CLI progress log lines."""

    return time.strftime("[%H:%M:%S] ")


def _print_list(
    species: Optional[str],
    annotation: Optional[str],
    isoform_set: Optional[str],
) -> None:
    print("species\tisoform_set\tannotation\tversion\tdefault\tinstalled_path\tsource_gene")
    for entry in iter_resources():
        if species is not None and entry.species != species:
            continue
        if annotation is not None and entry.annotation != annotation:
            continue
        if isoform_set is not None and entry.isoform_set != isoform_set:
            continue
        default = "yes" if entry.is_default else "no"
        print(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                entry.species,
                entry.isoform_set,
                entry.annotation,
                entry.version,
                default,
                entry.installed_relative_path,
                entry.source_relative_path,
            )
        )


if __name__ == "__main__":
    raise SystemExit(main())
