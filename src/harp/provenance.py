from __future__ import annotations

from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo
import hashlib
import json
import platform
import subprocess
from typing import Iterable


HARP_VERSION = "4.1"


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    path = Path(path)
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            h.update(chunk)
    return h.hexdigest()


def timestamp_now() -> datetime:
    return datetime.now().astimezone()


def timestamp_token(dt: datetime | None = None) -> str:
    dt = dt or timestamp_now()
    return dt.strftime("%Y%m%dT%H%M%S%z")


def panel_output_name(n_taxa: int, dt: datetime | None = None) -> str:
    if n_taxa < 1:
        raise ValueError("n_taxa must be >= 1")
    return f"HARP_v{HARP_VERSION}_panel_n{n_taxa}_{timestamp_token(dt)}"


def git_commit(repo: str | Path) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "HEAD"],
            text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        return None


def build_provenance(*, files: Iterable[str | Path] = (), repo: str | Path | None = None,
                     extra: dict | None = None) -> dict:
    now = timestamp_now()
    file_rows = []
    for value in files:
        p = Path(value).expanduser().resolve()
        if p.is_file():
            file_rows.append({
                "path": str(p),
                "size_bytes": p.stat().st_size,
                "sha256": sha256_file(p),
            })
    result = {
        "harp_version": HARP_VERSION,
        "run_timestamp": now.isoformat(),
        "timestamp_token": timestamp_token(now),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "git_commit": git_commit(repo) if repo else None,
        "input_files": file_rows,
    }
    if extra:
        result.update(extra)
    return result


def write_provenance(path: str | Path, provenance: dict) -> None:
    Path(path).write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
