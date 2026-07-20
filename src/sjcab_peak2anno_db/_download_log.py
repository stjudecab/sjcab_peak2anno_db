"""Helpers for recording downloaded source URLs."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from ._registry import PathLike, user_data_dir

DOWNLOAD_LOG_NAME = "download_urls.log"


def record_download_url(
    url: str,
    data_dir: Optional[PathLike] = None,
    destination: Optional[PathLike] = None,
) -> Path:
    """Append a downloaded URL record under the configured data directory."""

    root = user_data_dir(data_dir)
    root.mkdir(parents=True, exist_ok=True)
    log_path = root / DOWNLOAD_LOG_NAME

    fields = [datetime.now(timezone.utc).isoformat(), url]
    if destination is not None:
        fields.append(str(Path(destination).expanduser()))

    with log_path.open("a", encoding="utf-8") as handle:
        handle.write("{}\n".format("\t".join(fields)))
    return log_path
