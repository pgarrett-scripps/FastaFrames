"""
This module implements the core functions of filterframes.
"""

from dataclasses import dataclass, asdict
from io import TextIOWrapper, StringIO
from typing import Union, TextIO, List, Dict, Tuple, Generator, Iterable
from enum import Enum
import warnings

import pandas as pd

from .util import get_lines, convert_to_best_datatype

COLS = [
    "db",
    "unique_identifier",
    "entry_name",
    "protein_name",
    "organism_name",
    "organism_identifier",
    "gene_name",
    "protein_existence",
    "sequence_version",
    "protein_sequence",
]


class FastaFields(Enum):
    """
    An enum representing the FASTA fields.

    :ivar PROTEIN_NAME: Field for protein name.
    :ivar ORGANISM_NAME: Field for organism name.
    :ivar ORGANISM_ID: Field for organism ID.
    :ivar GENE_NAME: Field for gene name.
    :ivar PROTEIN_EXISTENCE: Field for protein existence.
    :ivar SEQUENCE_VERSION: Field for sequence version.
    """

    PROTEIN_NAME = "PN"
    ORGANISM_NAME = "OS"
    ORGANISM_ID = "OX"
    GENE_NAME = "GN"
    PROTEIN_EXISTENCE = "PE"
    SEQUENCE_VERSION = "SV"


@dataclass
class FastaEntry:
    """
    A dataclass representing a FASTA entry.

    :param db: Database source. Default is empty string.
    :type db: str
    :param unique_identifier: Unique identifier for the entry. Default is empty string.
    :type unique_identifier: str
    :param entry_name: Name of the entry. Default is empty string.
    :type entry_name: str
    :param protein_name: Name of the protein. Default is None.
    :type protein_name: str, optional
    :param organism_name: Name of the organism. Default is None.
    :type organism_name: str, optional
    :param organism_identifier: Identifier for the organism. Default is None.
    :type organism_identifier: str, optional
    :param gene_name: Name of the gene. Default is None.
    :type gene_name: str, optional
    :param protein_existence: Existence status of the protein. Default is None.
    :type protein_existence: str, optional
    :param sequence_version: Version of the sequence. Default is None.
    :type sequence_version: str, optional
    :param protein_sequence: Sequence of the protein. Default is empty string.
    :type protein_sequence: str
    """

    db: str = ""
    unique_identifier: str = ""
    entry_name: str = ""
    protein_name: str = None
    organism_name: str = None
    organism_identifier: str = None
    gene_name: str = None
    protein_existence: str = None
    sequence_version: str = None
    protein_sequence: str = ""

    def serialize(self) -> str:
        """
        Serializes the FastaEntry object to a FASTA string.

        :return: Serialized FASTA string.
        :rtype: str
        """

        optional_fields = [
            (f" {FastaFields.PROTEIN_NAME.value}=", self.protein_name),
            (f" {FastaFields.ORGANISM_NAME.value}=", self.organism_name),
            (f" {FastaFields.ORGANISM_ID.value}=", self.organism_identifier),
            (f" {FastaFields.GENE_NAME.value}=", self.gene_name),
            (f" {FastaFields.PROTEIN_EXISTENCE.value}=", self.protein_existence),
            (f" {FastaFields.SEQUENCE_VERSION.value}=", self.sequence_version),
        ]

        fasta_header = f">{self.db}|{self.unique_identifier}|{self.entry_name}"
        fasta_header += "".join(
            f"{key}{value}" for key, value in optional_fields if value
        )

        return f"{fasta_header}\n{self.protein_sequence}\n"


def fasta_to_entries(
    data: Union[str, TextIOWrapper, StringIO, TextIO], skip_error: bool = False
) -> Generator[FastaEntry, None, None]:
    """
    Converts FASTA content to a list of FastaEntry objects.

    :param data: FASTA content or a file-like object containing FASTA content.
    :type data: Union[str, TextIOWrapper, StringIO, TextIO]
    :param skip_error: If True, skips invalid FASTA entries instead of raising an error.
    :type skip_error: bool

    :return: A generator that yields FastaEntry objects.
    :rtype: Generator[FastaEntry, None, None]
    """

    # Ensure data is iterable line by line
    lines = get_lines(data)

    current_entry = None

    for line in lines:
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):  # new protein
            if current_entry is not None:
                yield current_entry
                current_entry = None

            try:
                current_entry = _fasta_str_to_entry(line)
            except ValueError as e:
                if skip_error:
                    current_entry = None
                    continue
                raise e

        elif current_entry:
            current_entry.protein_sequence += line

    if current_entry:
        yield current_entry


def entries_to_df(entries: Iterable[FastaEntry]) -> pd.DataFrame:
    """
    Converts a list of FastaEntry objects to a pandas DataFrame.

    :param entries: List of FastaEntry objects.
    :type entries: List[FastaEntry]

    :return: A pandas DataFrame representing the FastaEntry objects.
    :rtype: pd.DataFrame
    """

    fasta_df = pd.DataFrame([asdict(entry) for entry in entries])

    for col_name in fasta_df:
        fasta_df[col_name] = convert_to_best_datatype(fasta_df[col_name])
    return fasta_df


def to_df(
    data: Union[str, TextIOWrapper, StringIO, TextIO, List[FastaEntry]],
    skip_error: bool = False,
) -> pd.DataFrame:
    """
    Converts a FASTA input or list of FastaEntry objects to a pandas DataFrame.

    :param data: FASTA content, a file-like object containing FASTA content, or a list of FastaEntry objects.
    :type data: Union[str, TextIOWrapper, StringIO, TextIO, List[FastaEntry]]
    :param skip_error: If True, skips invalid FASTA entries instead of raising an error.
    :type skip_error: bool

    :return: A pandas DataFrame representing the FASTA content or FastaEntry objects.
    :rtype: pd.DataFrame
    """

    if isinstance(data, list):
        return entries_to_df(data)

    return entries_to_df(fasta_to_entries(data, skip_error))


def df_to_entries(df: pd.DataFrame) -> List[FastaEntry]:
    """
    Converts a fasta dataframe to a list of FastaEntry objects.

    :param df: The fasta dataframe.
    :type df: pd.DataFrame

    :return: List of FastaEntry objects.
    :rtype: List[FastaEntry]
    """

    entries = [FastaEntry(**row.to_dict()) for _, row in df[COLS].iterrows()]
    return entries


def entries_to_fasta(
    entries: Iterable[FastaEntry], output_file: str = None
) -> Union[StringIO, None]:
    """
    Converts a list of FastaEntry objects to a StringIO object or file containing the fasta content.

    :param entries: The list containing FastaEntry objects.
    :type entries: Iterable[FastaEntry]
    :param output_file: The path to the output file. If None, the function will return a StringIO.
    :type output_file: str, optional

    :return: A StringIO object containing the fasta content or None if an output file is provided.
    :rtype: Union[StringIO, None]
    """

    fasta_string = StringIO()
    for entry in entries:
        fasta_string.write(entry.serialize())

    fasta_string.seek(0)

    if output_file is not None:
        with open(file=output_file, mode="w", encoding="UTF-8") as f:
            f.write(fasta_string.getvalue())
        return None

    return fasta_string


def to_fasta(
    data: Union[pd.DataFrame, Iterable[FastaEntry]], output_file: str = None
) -> Union[StringIO, None]:
    """
    Converts a fasta dataframe or list of FastaEntries to a StringIO object or file containing the fasta content.

    :param data: The fasta dataframe or a list of FastaEntry objects.
    :type data: Union[pd.DataFrame, Iterable[FastaEntry]]
    :param output_file: The path to the output file. If None, the function will return a StringIO.
    :type output_file: str, optional
    :return: A StringIO object containing the fasta content or None if an output file is provided.
    :rtype: Union[StringIO, None]
    """

    if isinstance(data, pd.DataFrame):
        return entries_to_fasta(df_to_entries(data), output_file)

    return entries_to_fasta(data, output_file)


def _extract_fasta_header_elements(entry_str: str) -> List[str]:
    """
    Extracts the elements from the header line of a fasta entry.

    :param entry_str: The header line of a fasta entry.
    :type entry_str: str
    :return: List of elements extracted from the header line.
    :rtype: List[str]
    """

    line_elements = entry_str.rstrip().replace(">", "").split(" ")
    return line_elements


def _extract_initial_info(line_elements: List[str]) -> Tuple[str, str, str]:
    """
    Extracts the initial information, such as database, unique identifier, and entry name from the list of elements.

    :param line_elements: List of elements extracted from the header line of a fasta entry.
    :type line_elements: List[str]
    :return: Tuple containing database, unique identifier, and entry name.
    :rtype: Tuple[str, str, str]
    """

    first_element_parts = line_elements[0].split("|")

    if len(first_element_parts) == 3:
        db = first_element_parts[0]
        unique_identifier = first_element_parts[1]
        entry_name = first_element_parts[2]
        return db, unique_identifier, entry_name

    if len(first_element_parts) >= 1:
        # write warning

        warnings.warn(
            f"Invalid fasta header format: {line_elements[0]}, using only the first part as unique identifier."
        )

        return None, line_elements[0], None

    raise ValueError(f"Invalid fasta header format: {line_elements[0]}")


def _process_line_elements(line_elements: List[str]) -> Dict[str, List[str]]:
    """
    Processes the list of line elements and groups them into a dictionary.

    The first element is ignored since this only contains the >db|UniqueIdentifier|EntryName string. After this, the
    next expected element is ProteinName which is not specified by the XX= notation. The remaining components will be
    specified by the XX= notation, and can be parsed accordingly.

    :param line_elements: List of elements extracted from the header line of a fasta entry.
    :type line_elements: List[str]
    :return: Dictionary containing grouped elements.
    :rtype: Dict[str, List[str]]
    """

    info = {}
    current_state = FastaFields.PROTEIN_NAME.value

    for elem in line_elements[1:]:

        if "=" in elem:
            current_state = elem[
                :2
            ]  # Assuming that the field keys are always two characters long
            elem = elem[3:]

        if current_state not in {field.value for field in FastaFields}:
            raise ValueError(
                f"Unexpected element: {current_state} encountered. Line: {line_elements}"
            )

        info.setdefault(current_state, []).append(elem)

    return info


def _fasta_str_to_entry(fasta_str: str) -> FastaEntry:
    """
    Extracts fasta information from the given fasta str and creates a FastaEntry object.

    :param fasta_str: The header line of a fasta entry.
    :type fasta_str: str
    :return: FastaEntry object containing the extracted information.
    :rtype: FastaEntry
    """

    line_elements = _extract_fasta_header_elements(fasta_str)
    db, unique_identifier, entry_name = _extract_initial_info(line_elements)
    info = _process_line_elements(line_elements)

    def _join_list_values(data: Dict[str, List[str]]) -> Dict[str, str]:
        return {k: " ".join(v) if v else None for k, v in data.items()}

    joined_info = _join_list_values(info)

    return FastaEntry(
        db=db,
        unique_identifier=unique_identifier,
        entry_name=entry_name,
        protein_name=joined_info.get("PN"),
        organism_name=joined_info.get("OS"),
        organism_identifier=joined_info.get("OX"),
        gene_name=joined_info.get("GN"),
        protein_existence=joined_info.get("PE"),
        sequence_version=joined_info.get("SV"),
    )
