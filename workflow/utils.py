from pathlib import Path
from urllib.parse import parse_qs, unquote, urlparse


def filename_from_url(url: str) -> str:
    """
    Extract a sensible filename from a URL, supporting DepMap API links that
    store the real name in the file_name query parameter.
    """
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    name = Path(parsed.path).name
    if (not name) or name == "download":
        file_names = query.get("file_name")
        if file_names:
            name = Path(unquote(file_names[0])).name
    name = unquote(name)
    if not name:
        raise ValueError(f"Unable to determine filename from URL: {url}")
    return name
