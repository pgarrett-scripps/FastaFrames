"""
This module implements the core functions of filterframes.
"""

from dataclasses import dataclass, asdict
from io import TextIOWrapper, StringIO
from typing import Union, TextIO, List, Dict, Tuple
from enum import Enum

import pandas as pd

from fastaframes.util import get_lines, convert_to_best_datatype


class FastaFields(Enum):
    """
    An enum representing the FASTA fields.
    """
    PROTEIN_NAME = 'PN'
    ORGANISM_NAME = 'OS'
    ORGANISM_ID = 'OX'
    GENE_NAME = 'GN'
    PROTEIN_EXISTENCE = 'PE'
    SEQUENCE_VERSION = 'SV'


@dataclass
class FastaEntry:
    """
    A dataclass representing a FASTA entry.
    """

    db: str = ''
    unique_identifier: str = ''
    entry_name: str = ''
    protein_name: str = None
    organism_name: str = None
    organism_identifier: str = None
    gene_name: str = None
    protein_existence: str = None
    sequence_version: str = None
    protein_sequence: str = ''

    def serialize(self) -> str:
        """
        Serializes the FastaEntry object to a FASTA string.
        """

        optional_fields = [
            (f' {FastaFields.PROTEIN_NAME.value}=', self.protein_name),
            (f' {FastaFields.ORGANISM_NAME.value}=', self.organism_name),
            (f' {FastaFields.ORGANISM_ID.value}=', self.organism_identifier),
            (f' {FastaFields.GENE_NAME.value}=', self.gene_name),
            (f' {FastaFields.PROTEIN_EXISTENCE.value}=', self.protein_existence),
            (f' {FastaFields.SEQUENCE_VERSION.value}=', self.sequence_version)
        ]

        fasta_header = f'>{self.db}|{self.unique_identifier}|{self.entry_name}'
        fasta_header += ''.join(f"{key}{value}" for key, value in optional_fields if value)

        return f"{fasta_header}\n{self.protein_sequence}\n"


def fasta_to_entries(data: Union[str, TextIOWrapper, StringIO, TextIO]) -> List[FastaEntry]:
    """
    Converts FASTA content to a list of FastaEntry objects.
    """

    lines = get_lines(data)

    entries = []
    for line in lines:
        if line == "":
            continue

        if line[0] == ">":  # new protein
            entries.append(_fasta_str_to_entry(line))
        else:
            entries[-1].protein_sequence += line.rstrip()

    return entries


def entries_to_df(entries: List[FastaEntry]) -> pd.DataFrame:
    """
    Converts a list of FastaEntry objects to a pandas DataFrame.
    """

    fasta_df = pd.DataFrame([asdict(entry) for entry in entries])
    for col_name in fasta_df:
        fasta_df[col_name] = convert_to_best_datatype(fasta_df[col_name])
    return fasta_df


def to_df(data: Union[str, TextIOWrapper, StringIO, TextIO, List[FastaEntry]]) -> pd.DataFrame:
    """
    Converts a FASTA input or list of FastaEntry objects to a pandas DataFrame.
    """

    if isinstance(data, list):
        return entries_to_df(data)

    return entries_to_df(fasta_to_entries(data))


def df_to_entries(df: pd.DataFrame) -> List[FastaEntry]:
    """
    Converts a fasta dataframe to a list of FastaEntry objects.
    """

    cols = ['db', 'unique_identifier', 'entry_name', 'protein_name', 'organism_name', 'organism_identifier',
            'gene_name', 'protein_existence', 'sequence_version', 'protein_sequence']

    entries = [FastaEntry(**row.to_dict()) for _, row in df[cols].iterrows()]
    return entries


def entries_to_fasta(entries: List[FastaEntry], output_file: str = None) -> Union[StringIO, None]:
    """
    Converts a list of FastaEntry objects to a StringIO object or file containing the fasta content.

    Args:
        entries (List[FastaEntry]): The list containing FastaEntry objects.
        output_file (str): The path to the output file, if None to_fasta will return a StringIO.

    Returns:
        StringIO: A StringIO object containing the fasta content.
    """

    fasta_string = StringIO()
    for entry in entries:
        fasta_string.write(entry.serialize())

    fasta_string.seek(0)

    if output_file is not None:
        with open(file=output_file, mode='w', encoding='UTF-8') as f:
            f.write(fasta_string.getvalue())
        return None

    return fasta_string


def to_fasta(data: Union[pd.DataFrame, List[FastaEntry]], output_file: str = None) -> Union[StringIO, None]:
    """
    Converts a fasta dataframe or list of FastaEntries to a StringIO object or file containing the fasta content.

    Args:
        data (pd.DataFrame): The fasta dataframe.
        output_file (str): The path to the output file, if None to_fasta will return a StringIO

    Returns:
        StringIO: A StringIO object containing the fasta content.
    """

    if isinstance(data, pd.DataFrame):
        return entries_to_fasta(df_to_entries(data), output_file)

    return entries_to_fasta(data, output_file)


def _extract_fasta_header_elements(entry_str: str) -> List[str]:
    """
    Extracts the elements from the header line of a fasta entry.
    """

    line_elements = entry_str.rstrip().replace('>', '').split(" ")
    return line_elements


def _extract_initial_info(line_elements: List[str]) -> Tuple[str, str, str]:
    """
    Extracts the initial information, such as database, unique identifier, and entry name from the list of elements.
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
    Processes the list of line elements and groups them into a dictionary. The first element is ignored since this
    only contains the >db|UniqueIdentifier|EntryName string. After this, the next expected element is ProteinName which
    is not specified by the XX= notation. The remaining components will be specified by the XX= notation, and can be
    parsed accordingly.

    Example Format:
    >db|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier GN=GeneName PE=ProteinExistence \
     SV=SequenceVersion

     Values are mapped as follows:
     {'PN':'ProteinName',
     'OS':'OrganismName',
     'OX':'OrganismIdentifier,
     'GN':'GeneName',,
     'PE':'ProteinExistence',,
     'SV':'SequenceVersion',}

    Args:
        line_elements (List[str]): A list of elements from the fasta header line.

    Returns:
        Dict[str, List[str]]: A dictionary with keys representing element categories
                              and values as lists of related elements.
    """

    info = {}
    current_state = FastaFields.PROTEIN_NAME.value  # Using Enum for readability
    for elem in line_elements[1:]:
        if '=' in elem:
            current_state = elem[:2]  # Assuming that the field keys are always two characters long
            elem = elem[3:]
        info.setdefault(current_state, []).append(elem)
    return info


def _fasta_str_to_entry(fasta_str: str) -> FastaEntry:
    """
    Extracts fasta information from the given fasta str and creates a FastaEntry object.
    """

    line_elements = _extract_fasta_header_elements(fasta_str)
    db, unique_identifier, entry_name = _extract_initial_info(line_elements)
    info = _process_line_elements(line_elements)

    def _join_list_values(data: Dict[str, List[str]]) -> Dict[str, str]:
        return {k: ' '.join(v) if v else None for k, v in data.items()}

    joined_info = _join_list_values(info)

    return FastaEntry(
        db=db,
        unique_identifier=unique_identifier,
        entry_name=entry_name,
        protein_name=joined_info.get('PN'),
        organism_name=joined_info.get('OS'),
        organism_identifier=joined_info.get('OX'),
        gene_name=joined_info.get('GN'),
        protein_existence=joined_info.get('PE'),
        sequence_version=joined_info.get('SV')
    )
