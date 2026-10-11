"""
Input normalization for fastaframes: turn any accepted FASTA source into something
fastatacular reads (a path or a file handle), or into a stream of lines.

Accepted sources: a ``str`` or :class:`os.PathLike` path (plain, gzip, bzip2 or xz),
raw FASTA text in a ``str``, a readable text or binary file-like object, or an
iterable of ``str``/``bytes`` lines.
"""

import bz2
import gzip
import io
import logging
import lzma
import os
from collections.abc import Generator, Iterable, Iterator
from io import StringIO
from pathlib import Path
from typing import IO, Any

from .exceptions import InvalidInputError

logger = logging.getLogger(__name__)


def _looks_like_text(s: str) -> bool:
    """A ``str`` is FASTA text, not a path, if it is blank, has a line break or starts with ``>``."""
    return not s.strip() or "\n" in s or "\r" in s or s.lstrip().startswith(">")


def as_source(data: Any) -> Path | IO[str] | IO[bytes]:
    """Return a path or open handle that :func:`fastatacular.read_fasta` accepts.

    - ``os.PathLike`` -> :class:`pathlib.Path`.
    - ``str`` that is FASTA text (blank, contains a newline, or starts with ``>``) ->
      :class:`io.StringIO`; any other ``str`` is a path.
    - a readable file-like object -> itself, rewound with ``seek(0)`` when seekable.
    - any other iterable of ``str``/``bytes`` lines -> :class:`io.StringIO` of the lines.

    :param data: The FASTA source.
    :return: A :class:`pathlib.Path` or a readable handle.
    :raises FileNotFoundError: If a path does not exist.
    :raises InvalidInputError: If ``data`` is none of the above (e.g. an ``int``).
    """
    if isinstance(data, str) and _looks_like_text(data):
        logger.debug("as_source: treating str input as raw FASTA text (%d chars)", len(data))
        return StringIO(data)
    if isinstance(data, (str, os.PathLike)):
        path = Path(data)
        if not path.exists():
            raise FileNotFoundError(f"FASTA file not found: {str(path)!r}")
        logger.debug("as_source: reading from file path %r", str(path))
        return path
    if hasattr(data, "read"):
        if getattr(data, "seekable", lambda: False)():
            data.seek(0)
        logger.debug("as_source: reading from file-like object %r", type(data).__name__)
        return data
    try:
        iterator = iter(data)
    except TypeError:
        logger.error("as_source: unsupported input type %r", type(data).__name__)
        raise InvalidInputError(data) from None
    logger.debug("as_source: reading from iterable of type %r", type(data).__name__)
    text = StringIO()
    for line in iterator:
        if isinstance(line, bytes):
            line = line.decode("UTF-8")
        text.write(line.rstrip("\r\n") + "\n")
    text.seek(0)
    return text


_MAGIC = ((b"\x1f\x8b", gzip.open), (b"BZh", bz2.open), (b"\xfd7zXZ\x00", lzma.open))


def _open_path(path: Path) -> IO[str]:
    """Open ``path`` as UTF-8 text, decompressing gzip/bzip2/xz by their magic bytes."""
    with open(path, "rb") as f:
        head = f.read(6)
    for magic, opener in _MAGIC:
        if head.startswith(magic):
            return opener(path, "rt", encoding="utf-8-sig")
    return open(path, encoding="utf-8-sig")


def get_lines(file_input: Any) -> Generator[str, None, None]:
    """Yield the lines of a FASTA source without their line endings.

    Accepts the same inputs as :func:`as_source`; bytes are decoded as UTF-8.

    :param file_input: The FASTA source.
    :return: A generator of lines.
    :raises FileNotFoundError: If a path does not exist.
    :raises InvalidInputError: If ``file_input`` is not a supported source.
    """
    source = as_source(file_input)
    handle = _open_path(source) if isinstance(source, Path) else source
    try:
        for line in handle:
            if isinstance(line, bytes):
                line = line.decode("UTF-8")
            yield line.rstrip("\r\n")
    finally:
        if handle is not source:
            handle.close()


class LineReader(io.TextIOBase):
    """A read-only text handle over an iterable of lines (each ending in ``\\n``).

    Lets :class:`fastatacular.FastaReader` stream lines that were cleaned up first,
    without reading the whole source into memory.
    """

    def __init__(self, lines: Iterable[str]) -> None:
        super().__init__()
        self._lines: Iterator[str] = iter(lines)

    def readable(self) -> bool:
        return True

    def __iter__(self) -> "LineReader":
        return self

    def __next__(self) -> str:
        return next(self._lines)

    def readline(self, size: int | None = -1, /) -> str:
        return next(self._lines, "")
