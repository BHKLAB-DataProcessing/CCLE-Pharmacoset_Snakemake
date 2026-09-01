#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path
from urllib.parse import urlparse


HTML_MARKERS = (
    b"<!DOCTYPE html",
    b"<html",
    b"DepMap \xe2\x80\x94 Verification",
    b"Quick check before you enter",
    b"Cookies must be enabled to access this site",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download or copy a resource with basic content validation."
    )
    parser.add_argument("--source", required=True, help="HTTP(S) URL or local path.")
    parser.add_argument("--output", required=True, help="Destination file path.")
    parser.add_argument(
        "--decompress-gzip",
        action="store_true",
        help="Validate the source as gzip and write decompressed output.",
    )
    parser.add_argument(
        "--expect-gzip",
        action="store_true",
        help="Validate the downloaded/copied output as gzip.",
    )
    parser.add_argument(
        "--expect-zip",
        action="store_true",
        help="Validate the downloaded/copied output as zip.",
    )
    return parser.parse_args()


def is_http_url(source: str) -> bool:
    return urlparse(source).scheme in {"http", "https"}


def run_curl(source: str, destination: Path) -> None:
    command = [
        "curl",
        "--fail",
        "--location",
        "--retry",
        "3",
        "--retry-delay",
        "2",
        "--connect-timeout",
        "30",
        "--output",
        str(destination),
        source,
    ]
    subprocess.run(command, check=True)


def copy_local(source: str, destination: Path) -> None:
    source_path = Path(source)
    if not source_path.exists():
        raise FileNotFoundError(f"Local source does not exist: {source}")
    shutil.copyfile(source_path, destination)


def reject_html(path: Path, source: str) -> None:
    head = path.read_bytes()[:4096]
    head_lower = head.lower()
    if any(marker.lower() in head_lower for marker in HTML_MARKERS):
        raise RuntimeError(
            "Downloaded content looks like an HTML verification/error page, "
            f"not raw data: {source}"
        )


def validate_gzip(path: Path, source: str) -> None:
    try:
        with gzip.open(path, "rb") as handle:
            handle.read(1)
    except OSError as exc:
        raise RuntimeError(f"Expected gzip content from {source}") from exc


def validate_zip(path: Path, source: str) -> None:
    if not zipfile.is_zipfile(path):
        raise RuntimeError(f"Expected zip content from {source}")


def write_decompressed_gzip(source_path: Path, output_path: Path) -> None:
    with gzip.open(source_path, "rb") as source, output_path.open("wb") as output:
        shutil.copyfileobj(source, output)


def main() -> int:
    args = parse_args()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=".download-", dir=output.parent) as tmp_dir:
        tmp_dir_path = Path(tmp_dir)
        downloaded = tmp_dir_path / "source"
        staged_output = tmp_dir_path / "output"

        print(f"[download] Fetching {args.source}", flush=True)
        if is_http_url(args.source):
            run_curl(args.source, downloaded)
        else:
            copy_local(args.source, downloaded)

        reject_html(downloaded, args.source)

        if args.decompress_gzip:
            validate_gzip(downloaded, args.source)
            write_decompressed_gzip(downloaded, staged_output)
        else:
            if args.expect_gzip:
                validate_gzip(downloaded, args.source)
            if args.expect_zip:
                validate_zip(downloaded, args.source)
            shutil.move(downloaded, staged_output)

        os.replace(staged_output, output)

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[download] ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
