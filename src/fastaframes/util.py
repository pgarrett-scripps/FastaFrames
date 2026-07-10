"""
Utility helpers for fastaframes: normalizing arbitrary inputs into lines.

Logs under the ``"fastaframes.util"`` logger at ``DEBUG`` so callers can trace
which input branch :func:`get_lines` took. See the package docstring for how to
enable logging.
"""

import logging
import os
from collections.abc import Generator
from io import StringIO, TextIOWrapper
from typing import TextIO

from .exceptions import InvalidInputError

logger = logging.getLogger(__name__)


def read_lines_from_text_io(file_input: TextIO) -> Generator[str, None, None]:
    """Yield newline-stripped lines from a readable file-like object.

    Handles both text and binary streams: ``bytes`` lines (e.g. from a stream
    opened in binary mode) are decoded as UTF-8. The stream is rewound with
    ``seek(0)`` first.
    """
    file_input.seek(0)
    for line in file_input:
        if isinstance(line, bytes):
            line = line.decode("UTF-8")
        yield line.rstrip("\n")


def read_lines_from_file(file_path: str) -> Generator[str, None, None]:
    """Yield newline-stripped lines from a UTF-8 encoded file on disk."""
    with open(file=file_path, encoding="UTF-8") as file:
        for line in file:
            yield line.rstrip("\n")


def read_lines_from_string(s: str) -> Generator[str, None, None]:
    """Yield lines from an in-memory string, splitting on ``"\\n"``."""
    for line in s.split("\n"):
        yield line.rstrip("\n")


def get_lines(
    file_input: str | TextIOWrapper | StringIO | TextIO,
) -> Generator[str, None, None]:
    """
    Normalize a FASTA source into a stream of newline-stripped lines.

    Dispatches on the input type:

    - ``str`` that is an existing path -> read the file.
    - ``str`` otherwise -> treat as raw FASTA text and split into lines.
    - a text file-like object (``TextIOWrapper``/``StringIO``) -> iterate it.
    - any other iterable of ``str``/``bytes`` lines -> iterate, decoding bytes.

    :param file_input: The FASTA source: a file path, a FASTA string, a readable
        file-like object, or an iterable of ``str``/``bytes`` lines.
    :return: A generator that yields each line without its trailing newline.
    :rtype: Generator[str, None, None]
    :raises InvalidInputError: If ``file_input`` is not a str, not a supported
        file-like object, and not iterable (e.g. an ``int``).
    """
    if isinstance(file_input, str):
        if os.path.exists(file_input):
            logger.debug("get_lines: reading from file path %r", file_input)
            yield from read_lines_from_file(file_input)
        else:
            logger.debug("get_lines: treating str input as raw FASTA text (%d chars)", len(file_input))
            yield from read_lines_from_string(file_input)
    elif isinstance(file_input, (TextIOWrapper, TextIO, StringIO)):
        logger.debug("get_lines: reading from text file-like object %r", type(file_input).__name__)
        yield from read_lines_from_text_io(file_input)
    else:
        try:
            iterator = iter(file_input)
        except TypeError:
            logger.error("get_lines: unsupported input type %r", type(file_input).__name__)
            raise InvalidInputError(file_input) from None
        logger.debug("get_lines: reading from iterable of type %r", type(file_input).__name__)
        for line in iterator:
            if isinstance(line, bytes):
                yield line.decode("UTF-8").rstrip("\n")
            else:
                yield line.rstrip("\n")
