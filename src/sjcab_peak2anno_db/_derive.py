"""Generate derived annotation BED files from bundled GeneBED files."""

from __future__ import annotations

from pathlib import Path
from typing import List, Tuple, Union

PathLike = Union[str, Path]


def write_tss(gene_bed: PathLike, output_bed: PathLike) -> Path:
    """Write 1 bp TSS intervals from a GeneBED file."""

    return _write_site(gene_bed, output_bed, site="tss")


def write_tes(gene_bed: PathLike, output_bed: PathLike) -> Path:
    """Write 1 bp TES intervals from a GeneBED file."""

    return _write_site(gene_bed, output_bed, site="tes")


def write_deduplong(
    gene_bed: PathLike,
    output_bed: PathLike,
    gene_key: str = "symbol",
) -> Path:
    """Keep the longest isoform per gene.

    By default, isoforms are grouped by gene symbol from column 4. Set
    ``gene_key="ensid"`` to group by the Ensembl/GENCODE gene ID from column 7
    with version suffix removed. Isoform length is read from column 5; if
    column 5 is not numeric, the BED interval length is used as a fallback.
    """

    key_mode = _normalize_gene_key(gene_key)
    output = Path(output_bed)
    output.parent.mkdir(parents=True, exist_ok=True)

    selected = {}
    with Path(gene_bed).open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = _parse_bed_fields(line, line_number)
            gene_name = _gene_key(fields, key_mode)
            length = _length(fields)
            previous = selected.get(gene_name)
            if previous is None or length > previous[0]:
                selected[gene_name] = (length, line)

    with output.open("w", encoding="utf-8") as handle:
        for _length_value, line in selected.values():
            handle.write(line)

    return output


def write_derived(
    gene_bed: PathLike,
    annotation: str,
    output_bed: PathLike,
    gene_key: str = "symbol",
) -> Path:
    """Write one derived annotation type from a GeneBED file."""

    annotation = annotation.lower()
    if annotation == "tss":
        return write_tss(gene_bed, output_bed)
    if annotation == "tes":
        return write_tes(gene_bed, output_bed)
    if annotation in {"deduplong", "dedup-long", "dedup_long"}:
        return write_deduplong(gene_bed, output_bed, gene_key=gene_key)
    raise ValueError("Unsupported derived annotation: {!r}".format(annotation))


def _write_site(gene_bed: PathLike, output_bed: PathLike, site: str) -> Path:
    output = Path(output_bed)
    output.parent.mkdir(parents=True, exist_ok=True)

    with Path(gene_bed).open("r", encoding="utf-8") as source, output.open(
        "w", encoding="utf-8"
    ) as target:
        for line_number, line in enumerate(source, start=1):
            if not line.strip() or line.startswith("#"):
                target.write(line)
                continue
            fields = _parse_bed_fields(line, line_number)
            fields[1], fields[2] = _site_interval(fields, site)
            target.write("\t".join(fields) + "\n")

    return output


def _parse_bed_fields(line: str, line_number: int) -> List[str]:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 6:
        raise ValueError(
            "Expected at least 6 BED columns at line {}; found {}.".format(
                line_number, len(fields)
            )
        )
    try:
        start = int(fields[1])
        end = int(fields[2])
    except ValueError as exc:
        raise ValueError(
            "Expected integer BED start/end at line {}.".format(line_number)
        ) from exc
    if start < 0 or end < start:
        raise ValueError(
            "Invalid BED interval at line {}: {}-{}.".format(
                line_number, start, end
            )
        )
    if fields[5] not in {"+", "-", "."}:
        raise ValueError(
            "Expected strand '+', '-', or '.' at line {}; found {!r}.".format(
                line_number, fields[5]
            )
        )
    return fields


def _site_interval(fields: List[str], site: str) -> Tuple[str, str]:
    start = int(fields[1])
    end = int(fields[2])
    strand = fields[5]
    if site == "tss":
        point = end - 1 if strand == "-" else start
    elif site == "tes":
        point = start if strand == "-" else end - 1
    else:
        raise ValueError("Unsupported site: {!r}".format(site))
    point = max(point, 0)
    return str(point), str(point + 1)


def _length(fields: List[str]) -> float:
    try:
        return float(fields[4])
    except ValueError:
        return float(int(fields[2]) - int(fields[1]))


def _gene_key(fields: List[str], gene_key: str) -> str:
    if gene_key == "ensid" and len(fields) > 6 and fields[6]:
        return _strip_version(fields[6])
    return fields[3]


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


def _strip_version(identifier: str) -> str:
    return identifier.split(".", 1)[0] if identifier else ""
