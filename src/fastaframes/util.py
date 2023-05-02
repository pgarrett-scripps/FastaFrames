"""
This module provides some simple utility functions for filterframes
"""

import os
from io import TextIOWrapper, StringIO
from typing import Union, List, TextIO, Any


def get_lines(file_input: Union[str, TextIOWrapper, StringIO, TextIO]) -> List[str]:
    """
    This function reads lines from a given input, which can be a file path, a string containing lines,
    a TextIOWrapper, or a StringIO object.

    Args:
        file_input (Union[str, TextIOWrapper, StringIO]): The input source.

    Returns:
        list: A list of lines from the input source.

    Raises:
        ValueError: If the input type is not supported.
    """
    if isinstance(file_input, str):
        if os.path.exists(file_input):
            with open(file=file_input, mode='r', encoding='UTF-8') as file:
                lines = file.read().split('\n')
        else:
            lines = file_input.split('\n')
    elif isinstance(file_input, (TextIOWrapper, TextIO)):
        lines = file_input.read().split('\n')

    elif isinstance(file_input, StringIO):
        lines = file_input.getvalue().split('\n')
    else:
        raise ValueError(f'Unsupported input type: {type(file_input)}!')

    return lines


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
