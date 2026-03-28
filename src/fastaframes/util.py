"""
This module provides utility functions for fastaframes.
"""

import os
from collections.abc import Generator
from io import StringIO, TextIOWrapper
from typing import TextIO


def read_lines_from_text_io(file_input: TextIO) -> Generator[str, None, None]:
    """Read lines from a TextIO object."""
    file_input.seek(0)
    for line in file_input:
        yield line.rstrip("\n")


def read_lines_from_file(file_path: str) -> Generator[str, None, None]:
    """Read lines from a file."""
    with open(file=file_path, encoding="UTF-8") as file:
        for line in file:
            yield line.rstrip("\n")


def read_lines_from_string(s: str) -> Generator[str, None, None]:
    """Read lines from a string."""
    for line in s.split("\n"):
        yield line.rstrip("\n")


def get_lines(
    file_input: str | TextIOWrapper | StringIO | TextIO,
) -> Generator[str, None, None]:
    """
    Retrieve lines from a file or string input.

    This function reads lines from a given input, which can be a file path, a string containing lines,
    a TextIOWrapper, or a StringIO object.

    Args:
        file_input: The input source.

    Returns:
        A generator that yields lines from the input source.

    Raises:
        ValueError: If the input type is not supported.
    """
    if isinstance(file_input, str):
        if os.path.exists(file_input):
            yield from read_lines_from_file(file_input)
        else:
            yield from read_lines_from_string(file_input)
    elif isinstance(file_input, (TextIOWrapper, TextIO, StringIO)):
        yield from read_lines_from_text_io(file_input)
    else:
        for line in file_input:
            if isinstance(line, bytes):
                yield line.decode("UTF-8").rstrip("\n")
            else:
                yield line.rstrip("\n")
