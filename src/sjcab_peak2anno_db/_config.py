"""User configuration for the sjcab peak-to-annotation database."""

from __future__ import annotations

import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

CONFIG_ENV_PREFIX = "SJCAB_PEAK2ANNO_DB"
DEFAULT_INSTALL_COMPONENTS = ("anno-bed", "anno-feature", "blacklists", "cgi")
DEFAULT_FEATURE_SPECS = (
    ("hg38", "v31"),
    ("hg19", "v31lift37"),
    ("mm10", "vM22"),
    ("mm39", "vM39"),
)
DEFAULT_STALE_DAYS = 90
DEFAULT_CLEAN_CACHE_DAYS = 90


@dataclass(frozen=True)
class UserConfig:
    """Resolved configuration values."""

    db_path: Optional[str] = None
    install_components: Tuple[str, ...] = DEFAULT_INSTALL_COMPONENTS
    feature_specs: Tuple[Tuple[str, str], ...] = DEFAULT_FEATURE_SPECS
    stale_days: int = DEFAULT_STALE_DAYS
    clean_cache_days: int = DEFAULT_CLEAN_CACHE_DAYS
    sizes_clean: bool = True
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
    specs_value = (
        values.get("install_species_versions")
        or values.get("default_install_species_versions")
        or values.get("species_versions")
    )
    if specs_value is None:
        species_value = values.get("install_species") or values.get("default_install_species")
        version_value = values.get("install_versions") or values.get("default_install_versions")
        if species_value and version_value:
            specs_value = _pair_lists(species_value, version_value)
    components = _parse_components(components_value)
    specs = _parse_specs(specs_value)
    stale_text = (
        values.get("species_txt_stale_days")
        or values.get("version_species_stale_days")
        or values.get("version_species_txt_stale_days")
        or values.get("ensembl_stale_days")
        or values.get("stale_days")
    )
    stale_days = _parse_stale_days(stale_text)
    clean_cache_text = values.get("cleancache") or values.get("clean_cache")
    clean_cache_days = _parse_clean_cache_days(clean_cache_text)
    sizes_clean = _parse_bool(
        values.get("sizesclean") or values.get("sizes_clean"), default=True
    )
    return UserConfig(
        db_path=db_path,
        install_components=components,
        feature_specs=specs,
        stale_days=stale_days,
        clean_cache_days=clean_cache_days,
        sizes_clean=sizes_clean,
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
    paths = [path for path in (home, xdg) if path.is_file()]
    explicit = os.environ.get("SJCAB_PEAK2ANNO_CONFIG")
    if explicit:
        explicit_path = Path(explicit).expanduser()
        if explicit_path.is_file():
            paths.append(explicit_path)
    # An explicitly selected file is the highest-precedence RC file.
    return tuple(paths)


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
    aliases = {"bed": "anno-bed", "feature": "anno-feature"}
    components = []
    for item in value.replace(";", ",").split(","):
        component = item.strip().lower()
        if component:
            components.append(aliases.get(component, component))
    return tuple(components) or DEFAULT_INSTALL_COMPONENTS


def _parse_specs(value: Optional[str]) -> Tuple[Tuple[str, str], ...]:
    if not value:
        return DEFAULT_FEATURE_SPECS
    specs = []
    for item in value.replace(";", ",").split(","):
        item = item.strip()
        if not item:
            continue
        if ":" in item:
            species, version = item.split(":", 1)
        elif "=" in item:
            species, version = item.split("=", 1)
        else:
            raise ValueError(
                "Install species/version must use species:version: {!r}".format(item)
            )
        if species.strip() and version.strip():
            specs.append((species.strip(), version.strip()))
    return tuple(specs) or DEFAULT_FEATURE_SPECS


def _pair_lists(species_value: str, version_value: str) -> str:
    species = _split_list(species_value)
    versions = _split_list(version_value)
    if len(species) != len(versions):
        raise ValueError("Install species and version lists must have equal length.")
    return ",".join("{}:{}".format(item, version) for item, version in zip(species, versions))


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
