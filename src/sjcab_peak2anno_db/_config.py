"""User configuration for the sjcab peak-to-annotation database."""

from __future__ import annotations

import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

CONFIG_ENV_PREFIX = "SJCAB_PEAK2ANNO_DB"
DEFAULT_INSTALL_COMPONENTS = ("genebed", "feature", "blacklists", "cgi")
DEFAULT_FEATURE_SPECS = (
    ("hg38", "v31"),
    ("hg19", "v31lift37"),
    ("mm10", "vM22"),
    ("mm39", "vM39"),
)
DEFAULT_STALE_DAYS = 90
DEFAULT_CLEAN_CACHE_DAYS = 90
DEFAULT_PEAK_TXT_DELIMITER = ":-*/^;_%$,"
DEFAULT_BED_SCORE_COLUMN = 5
DEFAULT_TXT_SCORE_COLUMN = 2
_RC_DEFAULT_LINES = (
    "#SJCAB_PEAK2ANNO_DB_PATH=~/.sjcab_peak2anno_db",
    "#SJCAB_PEAK2ANNO_DB_INSTALL_OPTIONS=genebed,feature,blacklists,cgi",
    "#SJCAB_PEAK2ANNO_DB_INSTALL_SPECIES=hg38,hg19,mm10,mm39",
    "#SJCAB_PEAK2ANNO_DB_INSTALL_VERSIONS=v31,v31lift37,vM22,vM39",
    "#SJCAB_PEAK2ANNO_DB_VERSION_STALE_DAYS=90",
    "#SJCAB_PEAK2ANNO_DB_SIZESCLEAN=1",
    "#SJCAB_PEAK2ANNO_DB_CLEANCACHE=90",
    "#SJCAB_PEAK2ANNO_DB_PEAK_TXT_DELIMITER=:-*/^;_%$,",
    "#SJCAB_PEAK2ANNO_DB_BED_SCORE_COLUMN=5",
    "#SJCAB_PEAK2ANNO_DB_TXT_SCORE_COLUMN=2",
)
_RC_VARIABLE_NAMES = frozenset(
    line[1:].split("=", 1)[0].upper() for line in _RC_DEFAULT_LINES
)


@dataclass(frozen=True)
class UserConfig:
    """Resolved configuration values."""

    db_path: Optional[str] = None
    install_components: Tuple[str, ...] = DEFAULT_INSTALL_COMPONENTS
    feature_specs: Tuple[Tuple[str, str], ...] = DEFAULT_FEATURE_SPECS
    stale_days: int = DEFAULT_STALE_DAYS
    clean_cache_days: int = DEFAULT_CLEAN_CACHE_DAYS
    sizes_clean: bool = True
    peak_txt_delimiter: str = DEFAULT_PEAK_TXT_DELIMITER
    bed_score_column: int = DEFAULT_BED_SCORE_COLUMN
    txt_score_column: int = DEFAULT_TXT_SCORE_COLUMN
    components_configured: bool = False
    feature_specs_configured: bool = False


def load_config() -> UserConfig:
    """Load RC files, then apply ``SJCAB_PEAK2ANNO_DB_*`` overrides."""

    values: Dict[str, str] = {}
    for path in _config_paths():
        values.update(_read_rc(path))

    env_values = {
        key[len(CONFIG_ENV_PREFIX) + 1 :].lower(): value
        for key, value in os.environ.items()
        if key.startswith(CONFIG_ENV_PREFIX + "_") and value.strip()
    }
    values.update(env_values)

    db_path = values.get("db_path") or values.get("path")
    components_value = values.get("install_options") or values.get("default_install_options")
    species_value = values.get("install_species") or values.get("default_install_species")
    version_value = values.get("install_versions") or values.get("default_install_versions")
    specs_value = None
    if species_value and version_value:
        specs_value = _pair_lists(species_value, version_value)
    components = _parse_components(components_value)
    specs = _parse_specs(specs_value)
    stale_text = (
        values.get("version_stale_days")
        or values.get("ensembl_stale_days")
        or values.get("stale_days")
    )
    stale_days = _parse_stale_days(stale_text)
    clean_cache_text = values.get("cleancache") or values.get("clean_cache")
    clean_cache_days = _parse_clean_cache_days(clean_cache_text)
    sizes_clean = _parse_bool(
        values.get("sizesclean") or values.get("sizes_clean"), default=True
    )
    peak_txt_delimiter = (
        values.get("peak_txt_delimiter") or DEFAULT_PEAK_TXT_DELIMITER
    ).strip()
    bed_score_column = _parse_column(
        values.get("bed_score_column"), DEFAULT_BED_SCORE_COLUMN
    )
    txt_score_column = _parse_column(
        values.get("txt_score_column"), DEFAULT_TXT_SCORE_COLUMN
    )
    return UserConfig(
        db_path=db_path,
        install_components=components,
        feature_specs=specs,
        stale_days=stale_days,
        clean_cache_days=clean_cache_days,
        sizes_clean=sizes_clean,
        peak_txt_delimiter=peak_txt_delimiter,
        bed_score_column=bed_score_column,
        txt_score_column=txt_score_column,
        components_configured=components_value is not None,
        feature_specs_configured=specs_value is not None,
    )


def configured_db_path() -> Optional[str]:
    return load_config().db_path


def configured_stale_seconds() -> int:
    return load_config().stale_days * 86400


def configured_clean_cache_days() -> int:
    return load_config().clean_cache_days


def configured_sizes_clean() -> bool:
    return load_config().sizes_clean


def configured_peak_txt_delimiter() -> str:
    return load_config().peak_txt_delimiter


def configured_bed_score_column() -> int:
    return load_config().bed_score_column


def configured_txt_score_column() -> int:
    return load_config().txt_score_column


def effective_processes(requested: int) -> int:
    """Limit workers to the CPUs allocated to the current process."""

    if requested < 1:
        raise ValueError("processes must be at least 1")
    try:
        allocated = len(os.sched_getaffinity(0))
    except AttributeError:
        return requested
    return max(1, min(requested, allocated))


def clean_cache_files(
    cache_dir: Path,
    policy=None,
    current_path: Optional[Path] = None,
) -> None:
    """Remove cache files according to an age policy.

    ``None`` uses the configured age, a negative value removes the current
    file immediately, and ``False`` disables cleanup for that call.
    """

    if policy is False:
        return
    if policy is None:
        days = configured_clean_cache_days()
    elif isinstance(policy, bool):
        days = -1 if policy else None
    else:
        days = int(policy)
    if days is None:
        return
    root = Path(cache_dir)
    if not root.is_dir():
        return
    if days < 0:
        candidates = (Path(current_path),) if current_path is not None else tuple()
    else:
        cutoff = time.time() - days * 86400
        candidates = tuple(
            path for path in root.rglob("*") if path.is_file() and path.stat().st_mtime < cutoff
        )
    for path in candidates:
        if path.is_file():
            path.unlink()


def _config_paths() -> Tuple[Path, ...]:
    home = Path.home() / ".sjcab_peak2anno.rc"
    xdg_root = Path(
        os.environ.get("XDG_CONFIG_HOME", str(Path.home() / ".config"))
    ).expanduser()
    xdg = xdg_root / "sjcab_peak2anno" / ".sjcab_peak2anno.rc"
    _migrate_legacy_default_xdg_config(xdg)
    explicit = os.environ.get("SJCAB_PEAK2ANNO_CONFIG")
    paths = [path for path in (home, xdg) if path.is_file()]
    if explicit:
        explicit_path = Path(explicit).expanduser()
        if explicit_path.is_file():
            paths.append(explicit_path)
    if not paths:
        _create_default_xdg_config(xdg)
        if xdg.is_file():
            paths.append(xdg)
    else:
        for path in paths:
            _ensure_rc_variables(path)
    # An explicitly selected file is the highest-precedence RC file.
    return tuple(paths)


def _create_default_xdg_config(path: Path) -> None:
    """Create the default XDG RC file when no user RC file exists."""

    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            """# sjcab_peak2anno_db configuration.
# Uncomment or edit values as needed.
"""
            + "\n".join(_RC_DEFAULT_LINES)
            + "\n",
            encoding="utf-8",
        )
    except OSError:
        # Configuration should remain usable on read-only or restricted homes.
        return


def _migrate_legacy_default_xdg_config(path: Path) -> None:
    """Comment settings from the previously generated active template."""

    if not path.is_file():
        return
    defaults = {
        "SJCAB_PEAK2ANNO_DB_INSTALL_OPTIONS=genebed,feature,blacklists,cgi",
        "SJCAB_PEAK2ANNO_DB_INSTALL_SPECIES=hg38,hg19,mm10,mm39",
        "SJCAB_PEAK2ANNO_DB_INSTALL_VERSIONS=v31,v31lift37,vM22,vM39",
        "SJCAB_PEAK2ANNO_DB_VERSION_STALE_DAYS=90",
        "SJCAB_PEAK2ANNO_DB_SIZESCLEAN=1",
        "SJCAB_PEAK2ANNO_DB_CLEANCACHE=90",
        "SJCAB_PEAK2ANNO_DB_PEAK_TXT_DELIMITER=:-*/^;_%$,",
        "SJCAB_PEAK2ANNO_DB_BED_SCORE_COLUMN=5",
        "SJCAB_PEAK2ANNO_DB_TXT_SCORE_COLUMN=2",
    }
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
        active = {line.strip() for line in lines if line.strip() and not line.lstrip().startswith("#")}
        if lines and lines[0].strip() == "# sjcab_peak2anno_db configuration." and active <= defaults:
            _create_default_xdg_config(path)
    except OSError:
        return


def _ensure_rc_variables(path: Path) -> None:
    """Append missing commented configuration variables to an RC file."""

    try:
        text = path.read_text(encoding="utf-8")
        present = set()
        normalized_lines = []
        changed = False
        for raw_line in text.splitlines():
            line = raw_line.strip()
            if line.startswith("#"):
                comment = line[1:].lstrip()
                comment_key = comment.split("=", 1)[0].strip().upper()
                if comment_key in _RC_VARIABLE_NAMES:
                    normalized_line = "#" + comment
                    if normalized_line != raw_line:
                        changed = True
                    raw_line = normalized_line
                    line = comment
                else:
                    line = comment
            if "=" not in line:
                normalized_lines.append(raw_line)
                continue
            key = line.split("=", 1)[0].strip().upper()
            if key in _RC_VARIABLE_NAMES:
                present.add(key)
            normalized_lines.append(raw_line)
        missing = [
            line for line in _RC_DEFAULT_LINES
            if line[1:].split("=", 1)[0].upper() not in present
        ]
        updated = "\n".join(normalized_lines)
        if text.endswith("\n"):
            updated += "\n"
        if missing:
            suffix = "" if not updated or updated.endswith("\n") else "\n"
            updated += suffix + "\n".join(missing) + "\n"
        if not changed and updated == text:
            return
        path.write_text(updated, encoding="utf-8")
    except OSError:
        return


def _read_rc(path: Path) -> Dict[str, str]:
    values: Dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith(("#", ";", "[")):
            continue
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        prefix = CONFIG_ENV_PREFIX + "_"
        if key.upper().startswith(prefix):
            key = key[len(prefix) :]
        values[key.strip().lower().replace("-", "_")] = value.strip()
    return values


def _parse_components(value: Optional[str]) -> Tuple[str, ...]:
    if not value:
        return DEFAULT_INSTALL_COMPONENTS
    components = []
    for item in value.replace(";", ",").split(","):
        component = item.strip().lower()
        if component:
            components.append(component)
    return tuple(components) or DEFAULT_INSTALL_COMPONENTS


def parse_species_version_values(
    species_value: str,
    version_value: Optional[str] = None,
) -> Tuple[Tuple[str, Optional[str]], ...]:
    """Expand species/version values and pair or broadcast them.

    Values may be comma-separated or point to ``.lst``/``.list`` files. A
    two-column file, or a species value containing ``species:version``,
    supplies explicit pairs.
    """

    species_items, species_pairs = _choice_entries(species_value)
    version_items, version_pairs = _choice_entries(version_value)
    if species_pairs and (version_value or version_pairs):
        raise ValueError("Species:version pairs cannot be combined with a version list.")
    if version_pairs:
        if species_items:
            raise ValueError("A two-column version list cannot be combined with species.")
        return tuple(version_pairs)
    if species_pairs:
        return tuple(species_pairs)
    if not species_items:
        return tuple()
    if not version_items:
        return tuple((species, None) for species in species_items)
    if len(species_items) == 1:
        return tuple((species_items[0], version) for version in version_items)
    if len(version_items) == 1:
        return tuple((species, version_items[0]) for species in species_items)
    if len(species_items) == len(version_items):
        return tuple(zip(species_items, version_items))
    raise ValueError(
        "Species and version lists must have the same number of values, "
        "or one list must contain a single value."
    )


def _choice_entries(value: Optional[str]):
    if not value:
        return [], []
    path = Path(value).expanduser()
    if path.suffix.lower() in (".lst", ".list") and path.is_file():
        raw_items = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
    else:
        raw_items = [item.strip() for item in value.replace(";", ",").split(",")]
    items = []
    pairs = []
    for item in raw_items:
        if not item or item.startswith(("#", ";")):
            continue
        fields = item.split()
        if len(fields) >= 2:
            pairs.append((fields[0], fields[1]))
        elif ":" in item or "=" in item:
            separator = ":" if ":" in item else "="
            species, version = item.split(separator, 1)
            if species.strip() and version.strip():
                pairs.append((species.strip(), version.strip()))
        else:
            items.append(item)
    return items, pairs


def _parse_specs(value: Optional[str]) -> Tuple[Tuple[str, str], ...]:
    if not value:
        return DEFAULT_FEATURE_SPECS
    specs = parse_species_version_values(value)
    if any(version is None for _, version in specs):
        raise ValueError("Install species/version must use species:version.")
    return tuple(specs) or DEFAULT_FEATURE_SPECS


def _pair_lists(species_value: str, version_value: str) -> str:
    specs = parse_species_version_values(species_value, version_value)
    return ",".join("{}:{}".format(species, version) for species, version in specs)


def _split_list(value: str):
    return tuple(item.strip() for item in value.replace(";", ",").split(",") if item.strip())


def _parse_stale_days(value: Optional[str]) -> int:
    if not value:
        return DEFAULT_STALE_DAYS
    try:
        days = int(value)
    except ValueError as exc:
        raise ValueError("SJCAB_PEAK2ANNO_DB stale days must be an integer.") from exc
    if days < 0:
        raise ValueError("SJCAB_PEAK2ANNO_DB stale days cannot be negative.")
    return days


def _parse_clean_cache_days(value: Optional[str]) -> int:
    if not value:
        return DEFAULT_CLEAN_CACHE_DAYS
    try:
        return int(value)
    except ValueError as exc:
        raise ValueError(
            "SJCAB_PEAK2ANNO_DB_CLEANCACHE must be an integer number of days."
        ) from exc


def _parse_bool(value: Optional[str], default: bool = False) -> bool:
    if value is None:
        return default
    return value.strip().lower() in {"1", "true", "yes", "on"}


def _parse_column(value: Optional[str], default: int) -> int:
    if not value:
        return default
    try:
        column = int(value)
    except ValueError as exc:
        raise ValueError("Selector score columns must be positive integers.") from exc
    if column < 1:
        raise ValueError("Selector score columns must be positive integers.")
    return column
