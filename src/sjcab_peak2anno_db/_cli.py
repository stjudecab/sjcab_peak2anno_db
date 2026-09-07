"""Command-line interface for sjcab_peak2anno_db."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Optional, Sequence

from ._chromhmm import download_chromhmm
from ._dedup import dedup_gencode_bed, filter_gencode_bed
from ._external import CGI_SPECIES, download_cgi, install_blacklists, install_cgi
from ._gencode import download_and_convert_gencode_gtf
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
    path_parser.add_argument("-d", "--data-dir", help="Generated annotation directory.")
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
    install_parser.add_argument("-d", "--data-dir", help="Generated annotation directory.")
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
    update_parser.add_argument("-d", "--data-dir", help="Generated annotation directory.")

    blacklist_parser = subparsers.add_parser(
        "install-blacklists",
        help="Install blacklist BEDs into the configured user data directory.",
    )
    blacklist_parser.add_argument("-d", "--data-dir", help="Generated annotation directory.")
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
    install_cgi_parser.add_argument("-d", "--data-dir", help="Generated annotation directory.")
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
        help="Genome build to install, or all for default feature builds.",
    )
    feature_install_parser.add_argument(
        "version",
        nargs="?",
        help="GENCODE version to install, or all for default versions.",
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
        choices=_GENCODE_FEATURE_SPECIES,
        help="Genome build to install.",
    )
    feature_install_parser.add_argument(
        "-v",
        "--version",
        dest="version_option",
        help="GENCODE version to install.",
    )
    feature_install_parser.add_argument(
        "-d", "--data-dir", help="Generated annotation directory."
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
    )

    bed_install_parser = subparsers.add_parser(
        "install-bed",
        help="Install bundled GENCODE BEDs and derived TSS/TES files.",
    )
    bed_install_parser.add_argument("-d", "--data-dir", help="Generated annotation directory.")
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
    cgi_parser.add_argument("-d", "--data-dir", help="Generated annotation directory.")
    cgi_parser.add_argument(
        "-s",
        "--species",
        choices=CGI_SPECIES,
        nargs="+",
        help="Species to download. Defaults to all supported CGI species.",
    )
    cgi_parser.add_argument(
        "-n",
        "--no-overwrite",
        action="store_true",
        help="Do not rewrite existing CGI files.",
    )

    gencode_parser = subparsers.add_parser(
        "download-bed",
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
        help="GENCODE version, for example v31, v31lift37, or vM23.",
    )
    gencode_parser.add_argument(
        "-o",
        "--output-dir",
        default=".",
        help="Directory for generated BED files. GTFs use the database cache.",
    )
    gencode_parser.add_argument("-d", "--data-dir", help="Database/cache directory.")
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
        action="store_true",
        help="Delete the downloaded GTF from cachegtf after conversion.",
    )
    gencode_parser.add_argument(
        "-b",
        "--output-bed",
        help="Explicit single generated gene BED path.",
    )
    gencode_parser.add_argument(
        "-C",
        "--no-with-type",
        action="store_true",
        help="Write the compact 8-column BED instead of *.gene.bed.withtype.",
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

    _add_gencode_feature_parser(
        subparsers,
        "download-feature",
        "Download/convert GENCODE and write merged feature BED classes.",
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
                    args.version,
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
            )
            print(target)
            return 0
        if args.command == "update":
            target = update_data(data_dir=args.data_dir)
            print(target)
            return 0
        if args.command == "install-blacklists":
            target = install_blacklists(
                data_dir=args.data_dir,
                overwrite=not args.no_overwrite,
            )
            print(target)
            return 0
        if args.command == "install-cgi":
            target = install_cgi(
                data_dir=args.data_dir,
                overwrite=not args.no_overwrite,
            )
            print(target)
            return 0
        if args.command == "install-feature":
            feature_species, feature_version = _gencode_feature_install_scope(args)
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
                include_type=not args.compact_bed,
                tes_bp=args.tes_bp,
                progress=_stderr_progress,
            )
            print(target)
            return 0
        if args.command == "install-bed":
            target = install_gencode_beds(
                data_dir=args.data_dir,
                overwrite=args.overwrite,
            )
            print(target)
            return 0
        if args.command == "download-cgi":
            target = download_cgi(
                data_dir=args.data_dir,
                species=args.species,
                overwrite=not args.no_overwrite,
                progress=_stderr_progress,
            )
            print(target)
            return 0
        if args.command == "download-bed":
            target = download_and_convert_gencode_gtf(
                args.species,
                args.version,
                args.output_dir,
                gtf_path=args.gtf_path,
                gtf_url=args.url,
                output_bed=args.output_bed,
                include_type=not args.no_with_type,
                gene_types=args.gene_type,
                overwrite=not args.no_overwrite,
                progress=_stderr_progress,
            )
            print(target)
            return 0
        if args.command == "download-feature":
            specs = _gencode_feature_specs(args.species, args.version)
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
                    include_type=not args.compact_bed,
                    tes_bp=args.tes_bp,
                    overwrite=not args.no_overwrite,
                    progress=_stderr_progress,
                    cache_dir=args.data_dir,
                    ucsc_annotation=args.ucsc_source,
                    clean_cache=args.clean_cache,
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
                version=args.version,
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
            targets = download_chromhmm(
                data_dir=chromhmm_data_dir,
                model=args.model,
                genome=args.genome,
                ids=args.ids,
                tissue=args.tissue,
                cellline=args.cellline,
                overwrite=not args.no_overwrite,
                progress=_stderr_progress,
            )
            for eid in sorted(targets):
                print("{}\t{}".format(eid, targets[eid]))
            return 0
        if args.command in {"download-segway", "install-segway"}:
            segway_data_dir = (
                args.output_dir if args.command == "download-segway" else args.data_dir
            )
            requested_genome = args.genome
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
    parser.add_argument("species", help="Genome build, for example hg38 or hg19.")
    parser.add_argument(
        "version",
        help="GENCODE version, for example v31, v31lift37, or vM23.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default=".",
        help="Directory for generated gene and feature BED files.",
    )
    parser.add_argument("-d", "--data-dir", help="Database/cache directory.")
    _add_feature_generation_arguments(parser, include_output_dir=False)


def _add_dedup_filter_parser(
    subparsers: argparse._SubParsersAction,
    name: str,
    help_text: str,
) -> None:
    parser = subparsers.add_parser(name, help=help_text)
    parser.add_argument(
        "species",
        nargs="?",
        help="Genome build to resolve under the data directory when --gene-bed is absent.",
    )
    parser.add_argument(
        "version",
        nargs="?",
        default="default",
        help="GENCODE version to resolve. Defaults to the species default.",
    )
    parser.add_argument(
        "-b",
        "--gene-bed",
        help="Input all-isoform gene BED. Overrides species/version lookup.",
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
        "--data-dir",
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
) -> None:
    if include_output_dir:
        parser.add_argument(
            "-o",
            "--output-dir",
            default=".",
            help="Directory for generated gene and feature BED files.",
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
        action="store_true",
        help="Delete the downloaded GTF from cachegtf after conversion.",
    )
    parser.add_argument(
        "-b",
        "--gene-bed",
        help=(
            "Use an existing gene BED for TSS-flank files instead of converting "
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
    parser.add_argument(
        "-c",
        "--compact-bed",
        action="store_true",
        help="When converting GTF, write the compact 8-column gene BED.",
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
    parser = subparsers.add_parser(name, help=help_text)
    if name == "download-chromhmm":
        parser.add_argument(
            "-o",
            "--output-dir",
            default=".",
            help="Directory for downloaded ChromHMM files.",
        )
    else:
        parser.add_argument(
            "-d", "--data-dir", help="Generated annotation directory."
        )
    parser.add_argument(
        "-m",
        "--model",
        choices=("15", "18", "25"),
        default="18",
        help="Roadmap ChromHMM state model. Defaults to 18.",
    )
    parser.add_argument(
        "-G",
        "--genome",
        choices=("hg19", "hg38"),
        default="hg19",
        help="Genome build for BED files. Defaults to hg19; hg38 uses lifted-over BEDs.",
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
    parser = subparsers.add_parser(name, help=help_text)
    if name == "download-segway":
        parser.add_argument(
            "-o",
            "--output-dir",
            default=".",
            help="Directory for downloaded Segway files.",
        )
    else:
        parser.add_argument(
            "-d", "--data-dir", help="Generated annotation directory."
        )
    parser.add_argument(
        "-G",
        "--genome",
        default="hg19",
        help=(
            "Genome build for BED files. Segway source files are hg19; other "
            "builds require --yes-liftover to write a CrossMap script."
        ),
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
            "For non-hg19 --genome requests, download hg19 files and write a "
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
    if not species or not version:
        raise ValueError(
            "install-feature requires species and version, or --all."
        )
    return species, version


def _gencode_feature_specs(species: str, version: str) -> tuple:
    if species.lower() == "all":
        if version.lower() != "all":
            raise ValueError("download-feature all requires version all.")
        return DEFAULT_GENCODE_FEATURE_SPECS
    if version.lower() == "all":
        return ((species, _gencode_default_feature_version(species)),)
    return ((species, version),)


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
        print(message, end="", file=sys.stderr, flush=True)
        return
    if re.fullmatch(r"\d+\.\.", message):
        print(message, end="", file=sys.stderr, flush=True)
        return
    if message == "done":
        print(message, file=sys.stderr, flush=True)
        return
    print(message, file=sys.stderr, flush=True)


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
