"""Download helpers with optional progress callbacks."""

from __future__ import annotations

import shutil
import urllib.request
from pathlib import Path
from typing import Callable, Optional

ProgressCallback = Callable[[str], None]
USER_AGENT = "sjcab-peak2anno-db"


def report_progress(
    progress: Optional[ProgressCallback],
    message: str,
) -> None:
    """Send a progress message when a callback is configured."""

    if progress is not None:
        progress(message)


def download_file(
    url: str,
    destination: Path,
    timeout: int = 120,
    progress: Optional[ProgressCallback] = None,
) -> None:
    """Download a URL to a file, reporting progress every 5% when possible."""

    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = destination.with_name(destination.name + ".tmp")
    label = destination.name

    report_progress(progress, "download {}:".format(label))
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": USER_AGENT,
            "Accept": "application/octet-stream,*/*",
        },
    )

    try:
        with urllib.request.urlopen(request, timeout=timeout) as response, tmp_path.open(
            "wb"
        ) as output:
            size = _content_length(response)
            if size is None or size <= 0:
                shutil.copyfileobj(response, output)
            else:
                _copy_with_percent_progress(response, output, size, label, progress)
        tmp_path.replace(destination)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()
    report_progress(progress, "done")


def _content_length(response) -> Optional[int]:
    try:
        value = response.headers.get("Content-Length")
    except AttributeError:
        return None
    if not value:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def _copy_with_percent_progress(
    response,
    output,
    size: int,
    label: str,
    progress: Optional[ProgressCallback],
) -> None:
    downloaded = 0
    next_percent = 5
    while True:
        chunk = response.read(1024 * 1024)
        if not chunk:
            break
        output.write(chunk)
        downloaded += len(chunk)
        percent = min(100, int(downloaded * 100 / size))
        while percent >= next_percent and next_percent <= 100:
            report_progress(
                progress,
                "{}..".format(next_percent),
            )
            next_percent += 5
