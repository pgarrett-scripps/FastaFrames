"""
This module provides some simple utility functions for filterframes
"""

import os
from io import TextIOWrapper, StringIO
from typing import Union, List, TextIO, Any, Generator


def read_lines_from_text_io(file_input: TextIO) -> Generator[str, None, None]:
    """Read lines from a TextIO object."""
    file_input.seek(0)
    for line in file_input:
        yield line.rstrip('\n')


def read_lines_from_file(file_path: str) -> Generator[str, None, None]:
    """Read lines from a file."""
    with open(file=file_path, mode='r', encoding='UTF-8') as file:
        for line in file:
            yield line.rstrip('\n')


def read_lines_from_string(s: str) -> Generator[str, None, None]:
    """Read lines from a string."""
    for line in s.split('\n'):
        yield line.rstrip('\n')


def get_lines(file_input: Union[str, TextIOWrapper, StringIO, TextIO]) -> Generator[str, None, None]:
    """
    Retrieve lines from a file or string input.

    This function reads lines from a given input, which can be a file path, a string containing lines,
    a TextIOWrapper, or a StringIO object.

    Args:
        file_input (Union[str, TextIOWrapper, StringIO, TextIO]): The input source.

    Returns:
        generator: A generator that yields lines from the input source.

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
                yield line.decode('UTF-8').rstrip('\n')
            else:
                yield line.rstrip('\n')


def best_datatype_for_list(values: List[Any]) -> Union[type, None]:
    """
    Determines the best data type that can be applied to all values in the list.

    Args:
        values (List[Any]): A list of values.

    Returns:
        type: The best data type for the list, or None if no common type is found.
    """
    for datatype in [int, float, str]:
        try:
            _ = [datatype(value) for value in values]
            return datatype
        except (ValueError, TypeError):
            continue
    return None


def convert_to_best_datatype(values: List[Any]) -> List[Any]:
    """
    Tries to convert a list of values to the most specific datatype possible: int, float, or str, in that order.

    Args:
        values (List[Any]): A list of values to be converted.

    Returns:
        List[Any]: A list of converted values.

    Raises:
        ValueError: If unable to convert values to any datatype.
    """
    best_datatype = best_datatype_for_list(values)

    if best_datatype is None:
        raise ValueError("Unable to convert values to any datatype")

    return [best_datatype(value) for value in values]
