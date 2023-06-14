"""
This module provides some simple utility functions for filterframes
"""

import os
from io import TextIOWrapper, StringIO
from typing import Union, List, TextIO, Any, Generator

FILE_TYPES = Union[str, TextIOWrapper, StringIO, TextIO]


def _get_lines(file_input: FILE_TYPES) -> Generator[str, None, None]:
    """
    Retrieve lines from a file or string input.

    This function reads lines from a given input, which can be a file path, a string containing lines,
    a TextIOWrapper, or a StringIO object.

    Args:
        file_input (Union[str, TextIOWrapper, StringIO]): The input source.

    Returns:
        generator: A generator that yields lines from the input source.

    Raises:
        ValueError: If the input type is not supported.
    """
    if isinstance(file_input, str):  # File path or string
        if os.path.exists(file_input):
            with open(file=file_input, mode='r', encoding='UTF-8') as file:
                for line in file:
                    yield line.rstrip('\n')
        else:
            for line in file_input.split('\n'):
                yield line.rstrip('\n')
    elif isinstance(file_input, (TextIOWrapper, TextIO)):  # TextIOWrapper or StringIO
        file_input.seek(0)
        for line in file_input:
            yield line.rstrip('\n')
    elif isinstance(file_input, StringIO):  # StringIO
        file_input.seek(0)
        for line in file_input.readlines():
            yield line.rstrip('\n')
    else:
        try:
            for line in file_input:
                yield line.decode('UTF-8').rstrip('\n')
        except Exception as e:
            raise ValueError(f'Unsupported input type: {type(file_input)}!')


def convert_to_best_datatype(values: List[Any]):
    """
    This function tries to convert a list of values to either float, int, or str datatypes, in that order.

    Args:
        values (List[Any]): A list of values to be converted.

    Returns:
        list: A list of converted values.

    Raises:
        ValueError: If unable to convert values to any datatype.
    """

    for datatype in [float, int, str]:
        try:
            converted_values = [datatype(value) for value in values]
            return converted_values
        except (ValueError, TypeError):
            continue
    raise ValueError("Unable to convert values to any datatype")
