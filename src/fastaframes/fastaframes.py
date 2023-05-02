"""
This module provides functions for working with FASTA files.
"""

import os
from dataclasses import dataclass, asdict
from io import TextIOWrapper, StringIO
from typing import Union, TextIO, List, Dict, Tuple, Any

import pandas as pd


def _get_lines(file_input: Union[str, TextIOWrapper, StringIO, TextIO]) -> List[str]:
    """
    Retrieve lines from a file or string input.

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


def _convert_to_best_datatype(values: List[Any]):
    """
    Convert a list of values to the most suitable datatype.

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


@dataclass
class FastaEntry:
    """
        A data class representing a single entry in a FASTA file.
    """

    db: str = ''
    unique_identifier: str = ''
    entry_name: str = ''
    protein_name: Union[str, None] = None
    organism_name: Union[str, None] = None
    organism_identifier: Union[str, None] = None
    gene_name: Union[str, None] = None
    protein_existence: Union[str, None] = None
    sequence_version: Union[str, None] = None
    protein_sequence: str = ''

    def to_fasta(self) -> str:
        """
        Converts a FastaEntry object to a FASTA-formatted string.

        Returns:
            str: The FASTA-formatted string.

        """
        fasta_header = f'>{self.db}|{self.unique_identifier}|{self.entry_name}'
        if self.protein_name:
            fasta_header += f' PN={self.protein_name}'
        if self.organism_name:
            fasta_header += f' OS={self.organism_name}'
        if self.organism_identifier:
            fasta_header += f' OX={self.organism_identifier}'
        if self.gene_name:
            fasta_header += f' GN={self.gene_name}'
        if self.protein_existence:
            fasta_header += f' PE={self.protein_existence}'
        if self.sequence_version:
            fasta_header += f' SV={self.sequence_version}'

        return fasta_header + '\n' + self.protein_sequence + '\n'


def fasta_to_entries(file_input: Union[str, TextIOWrapper, StringIO, TextIO]) -> List[FastaEntry]:
    """
    Converts a FASTA file to a list of FastaEntry objects.

    Args:
        file_input (Union[str, TextIOWrapper, StringIO, TextIO]): A string or file object containing the FASTA data.

    Returns:
        List[FastaEntry]: A list of FastaEntry objects.

    """
    lines = _get_lines(file_input)
    entries = []
    for line in lines:
        if line == "":
            continue

        if line[0] == ">":  # new protein
            entries.append(_extract_fasta_info(line))
        else:
            entries[-1].protein_sequence += line.rstrip()

    return entries


def fasta_to_df(file_input: Union[str, TextIOWrapper, StringIO, TextIO]) -> pd.DataFrame:
    """
        Converts a FASTA file to a pandas DataFrame.

        Args:
            file_input (Union[str, TextIOWrapper, StringIO, TextIO]): A string or file object containing the FASTA data.

        Returns:
            pd.DataFrame: A pandas DataFrame containing the FASTA data.

    """
    entries = fasta_to_entries(file_input)
    fasta_df = pd.DataFrame([asdict(entry) for entry in entries])
    for col_name in fasta_df:
        fasta_df[col_name] = _convert_to_best_datatype(fasta_df[col_name])
    return fasta_df


def df_to_entries(fasta_df: pd.DataFrame) -> List[FastaEntry]:
    """
    Converts a fasta dataframe to a list of FastaEntry objects.

    Args:
        fasta_df (pd.DataFrame): The fasta dataframe.

    Returns:
        List[FastaEntry]: A list of FastaEntry objects.
    """
    entries = [FastaEntry(**row.to_dict()) for _, row in fasta_df.iterrows()]
    return entries


def to_fasta(fasta_input: Union[pd.DataFrame, List[FastaEntry]]) -> StringIO:
    """
    Converts a fasta dataframe to a StringIO object containing the fasta content.

    Args:
        fasta_input (pd.DataFrame): The fasta dataframe.

    Returns:
        StringIO: A StringIO object containing the fasta content.
    """
    fasta_string = StringIO()

    if isinstance(fasta_input, pd.DataFrame):
        entries = df_to_entries(fasta_input)
    else:
        entries = fasta_input

    for entry in entries:
        fasta_string.write(f'>{entry.db}|{entry.unique_identifier}|{entry.entry_name}')
        if entry.protein_name:
            fasta_string.write(f' PN={entry.protein_name}')
        if entry.organism_name:
            fasta_string.write(f' OS={entry.organism_name}')
        if entry.organism_identifier:
            fasta_string.write(f' OX={entry.organism_identifier}')
        if entry.gene_name:
            fasta_string.write(f' GN={entry.gene_name}')
        if entry.protein_existence:
            fasta_string.write(f' PE={entry.protein_existence}')
        if entry.sequence_version:
            fasta_string.write(f' SV={entry.sequence_version}')
        fasta_string.write('\n')
        fasta_string.write(entry.protein_sequence + '\n')

    fasta_string.seek(0)
    return fasta_string


def _extract_fasta_header_elements(fasta_entry: str) -> List[str]:
    """
    Extracts the elements from the header line of a fasta entry.

    Args:
        fasta_entry (str): The header line of a fasta entry.

    Returns:
        List[str]: A list of elements found in the header line.
    """
    line_elements = fasta_entry.rstrip().replace('>', '').split(" ")
    return line_elements


def _extract_initial_info(line_elements: List[str]) -> Tuple[str, str, str]:
    """
    Extracts the initial information, such as database, unique identifier,
    and entry name from the list of elements.

    Args:
        line_elements (List[str]): A list of elements from the fasta header line.

    Returns:
        Tuple[str, str, str]: A tuple containing the database, unique identifier, and entry name.

    Raises:
        ValueError: If the fasta entry format is invalid.
    """
    first_element_parts = line_elements[0].split('|')
    if len(first_element_parts) >= 3:
        db = first_element_parts[0]
        unique_identifier = first_element_parts[1]
        entry_name = first_element_parts[2]

        return db, unique_identifier, entry_name

    raise ValueError("Invalid fasta entry format")


def _process_line_elements(line_elements: List[str]) -> Dict[str, List[str]]:
    """
    Processes the list of line elements and groups them into a dictionary.

    Args:
        line_elements (List[str]): A list of elements from the fasta header line.

    Returns:
        Dict[str, List[str]]: A dictionary with keys representing element categories
                              and values as lists of related elements.
    """
    info = {}
    current_state = 'PN'
    for elem in line_elements[1:]:
        if '=' in elem:
            current_state = elem[:2]
            elem = elem[3:]
        info.setdefault(current_state, []).append(elem)
    return info


def _join_list_values(info: Dict[str, List[str]]) -> Dict[str, Union[str, None]]:
    """
    Joins the list values of the info dictionary into strings or None if empty.

    Args:
        info (Dict[str, List[str]]): A dictionary with keys representing element categories
                                     and values as lists of related elements.

    Returns:
        Dict[str, Union[str, None]]: A dictionary with the same keys as the input dictionary,
                                      but with values as strings or None if empty.
    """
    return {k: ' '.join(info[k]) if len(info[k]) > 0 else None for k in info}


def _extract_fasta_info(fasta_entry: str) -> FastaEntry:
    """
    Extracts fasta information from the given fasta entry and creates a FastaEntry object.

    Args:
        fasta_entry (str): The header line of a fasta entry.

    Returns:
        FastaEntry: An object containing the extracted fasta information.
    """
    line_elements = _extract_fasta_header_elements(fasta_entry)
    db, unique_identifier, entry_name = _extract_initial_info(line_elements)
    info = _process_line_elements(line_elements)
    info = _join_list_values(info)

    return FastaEntry(db=db,
                      unique_identifier=unique_identifier,
                      entry_name=entry_name,
                      protein_name=info.get('PN'),
                      organism_name=info.get('OS'),
                      organism_identifier=info.get('OX'),
                      gene_name=info.get('GN'),
                      protein_existence=info.get('PE'),
                      sequence_version=info.get('SV'))
