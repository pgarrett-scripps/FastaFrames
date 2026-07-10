"""
Core fastaframes functionality: the :class:`FastaEntry` record and the functions
that convert between FASTA text, ``FastaEntry`` objects, and pandas DataFrames.

Errors
------
Malformed FASTA headers raise :class:`fastaframes.exceptions.FastaFormatError`
(a subclass of ``ValueError``). Unrecognized ``KEY=value`` header fields are
preserved in :attr:`FastaEntry.additional_fields` rather than raising.

Logging
-------
Logs under the ``"fastaframes.fastaframes"`` logger: normal progress (entry
counts, output targets) at ``DEBUG`` and tolerated-but-odd input at ``WARNING``.
The package installs a ``NullHandler``, so nothing is emitted until logging is
configured.
"""

import logging
import warnings
from collections.abc import Generator, Iterable
from dataclasses import asdict, dataclass, field
from enum import Enum
from io import StringIO, TextIOWrapper
from typing import TextIO

import pandas as pd

from .exceptions import FastaFormatError
from .util import get_lines

logger = logging.getLogger(__name__)

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
    An enum representing the standard UniProt FASTA header fields.

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
    """A single FASTA record: the header fields plus the amino-acid sequence.

    Mirrors the UniProt-style header
    ``>db|unique_identifier|entry_name ProteinName OS=... OX=... GN=... PE=... SV=...``.
    Unrecognized ``KEY=value`` header fields are captured in
    :attr:`additional_fields` and re-emitted on serialization, so an entry
    round-trips through parse -> serialize without losing information.

    :param db: Source database (``"sp"`` for Swiss-Prot, ``"tr"`` for TrEMBL).
    :param unique_identifier: Primary accession (e.g. ``"P12345"``).
    :param entry_name: Entry name (e.g. ``"CP2D7_HUMAN"``).
    :param protein_name: Recommended protein name.
    :param organism_name: Scientific organism name (``OS=``).
    :param organism_identifier: NCBI taxonomy id (``OX=``).
    :param gene_name: Gene name (``GN=``).
    :param protein_existence: Protein-existence evidence level (``PE=``).
    :param sequence_version: Sequence version (``SV=``).
    :param protein_sequence: The amino-acid sequence.
    :param additional_fields: Non-standard ``KEY=value`` header fields, preserved
        for round-tripping. Defaults to an empty dict.
    """

    db: str | None = ""
    unique_identifier: str = ""
    entry_name: str | None = ""
    protein_name: str | None = None
    organism_name: str | None = None
    organism_identifier: str | None = None
    gene_name: str | None = None
    protein_existence: str | None = None
    sequence_version: str | None = None
    protein_sequence: str = ""
    additional_fields: dict[str, str] = field(default_factory=dict)

    @property
    def protein_id(self) -> str:
        """The ``db|unique_identifier|entry_name`` identifier.

        Empty components and the literal string ``"None"`` are dropped, so the
        result never contains a spurious ``"None"`` token.

        :return: The pipe-joined identifier, or ``""`` if no components are set.
        :rtype: str
        """
        elems = [self.db, self.unique_identifier, self.entry_name]
        elems = [e for e in elems if e is not None and e != "" and e != "None"]
        return "|".join(elems)

    @property
    def header(self) -> str:
        """The FASTA header line for this entry (leading ``>`` included).

        Builds ``>db|unique_identifier|entry_name`` followed by any set optional
        fields as `` KEY=value`` (the standard fields then
        :attr:`additional_fields`). Missing ``db``/``unique_identifier``/
        ``entry_name`` render as empty strings (e.g. ``>|P12345|``) rather than
        the literal text ``"None"``.

        :return: The FASTA header line (no trailing newline).
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
        optional_fields += [(f" {key}=", value) for key, value in self.additional_fields.items()]

        fasta_header = f">{self.db or ''}|{self.unique_identifier or ''}|{self.entry_name or ''}"
        fasta_header += "".join(f"{key}{value}" for key, value in optional_fields if value)
        return fasta_header

    def serialize(self, line_wrapping: bool = True, max_sequence_length: int | None = None) -> str:
        """Serialize this entry to a FASTA string: ``header`` + sequence + newline.

        The sequence is wrapped only when ``line_wrapping`` is True and
        ``max_sequence_length`` is a positive int; otherwise the sequence is
        written on a single line. The returned string ends with a trailing
        newline, so :func:`entries_to_fasta` can concatenate records directly.

        :param line_wrapping: Whether wrapping is allowed at all. Default True.
        :param max_sequence_length: Line width to wrap the sequence at. ``None``
            (default) or a non-positive value means "do not wrap".
        :return: The serialized FASTA record (with a trailing newline).
        :rtype: str
        """
        sequence = self.protein_sequence or ""
        if line_wrapping and sequence and max_sequence_length is not None and max_sequence_length > 0:
            sequence = "\n".join(
                sequence[i : i + max_sequence_length] for i in range(0, len(sequence), max_sequence_length)
            )
        return f"{self.header}\n{sequence}\n"

    def to_dict(self) -> dict[str, object]:
        """Return this entry as a dict of all fields plus the computed ``protein_id``.

        :return: Dictionary representation of the entry (includes ``protein_id``).
        :rtype: dict[str, object]
        """
        d: dict[str, object] = asdict(self)
        d["protein_id"] = self.protein_id
        return d


def fasta_to_entries(
    data: str | TextIOWrapper | StringIO | TextIO, skip_error: bool = False
) -> Generator[FastaEntry, None, None]:
    """
    Parse FASTA content and lazily yield one :class:`FastaEntry` per record.

    :param data: A file path, FASTA string, or file-like object (see
        :func:`fastaframes.util.get_lines` for accepted inputs).
    :param skip_error: If True, a header that fails to parse is skipped instead
        of raising; if False (default) the error propagates.
    :return: A generator that yields :class:`FastaEntry` objects.
    :rtype: Generator[FastaEntry, None, None]
    :raises FastaFormatError: On a malformed header when ``skip_error`` is False.
    """

    # Ensure data is iterable line by line
    lines = get_lines(data)

    current_entry = None
    yielded = 0
    skipped = 0

    for raw_line in lines:
        line = raw_line.strip()

        if not line:
            continue

        if line.startswith(">"):  # new protein
            if current_entry is not None:
                yield current_entry
                yielded += 1
                current_entry = None

            try:
                current_entry = _fasta_str_to_entry(line)
            except FastaFormatError as err:
                if skip_error:
                    skipped += 1
                    logger.warning("Skipping malformed FASTA header: %s", err)
                    current_entry = None
                    continue
                raise

        elif current_entry:
            current_entry.protein_sequence += line

    if current_entry:
        yield current_entry
        yielded += 1

    logger.debug("fasta_to_entries: yielded %d entries, skipped %d malformed headers", yielded, skipped)


def entries_to_df(entries: Iterable[FastaEntry]) -> pd.DataFrame:
    """
    Build a DataFrame from :class:`FastaEntry` objects (one row per entry).

    Columns are the :data:`COLS` fields plus ``additional_fields`` and the
    computed ``protein_id`` (see :meth:`FastaEntry.to_dict`).

    :param entries: The :class:`FastaEntry` objects to tabulate.
    :return: A pandas DataFrame representing the FastaEntry objects.
    :rtype: pd.DataFrame
    """

    fasta_df = pd.DataFrame([entry.to_dict() for entry in entries])
    return fasta_df


def to_df(
    data: str | TextIOWrapper | StringIO | TextIO | list[FastaEntry],
    skip_error: bool = False,
) -> pd.DataFrame:
    """
    Convert a FASTA source (or list of entries) into a pandas DataFrame.

    :param data: FASTA content, a file-like object, or a list of
        :class:`FastaEntry` objects.
    :param skip_error: If True, malformed FASTA records are skipped instead of
        raising. Ignored when ``data`` is already a list of entries.
    :return: A pandas DataFrame representing the FASTA content or entries.
    :rtype: pd.DataFrame
    :raises FastaFormatError: On malformed FASTA when ``skip_error`` is False.
    """

    if isinstance(data, list):
        return entries_to_df(data)  # ty: ignore[invalid-argument-type]

    return entries_to_df(fasta_to_entries(data, skip_error))


def df_to_entries(df: pd.DataFrame) -> list[FastaEntry]:
    """
    Convert a FASTA DataFrame into a list of :class:`FastaEntry` objects.

    Only the columns in :data:`COLS` are used; extra columns (such as
    ``protein_id``) are ignored, and any missing expected columns are left at
    their defaults with a :class:`UserWarning` naming them.

    :param df: The FASTA DataFrame (e.g. produced by :func:`to_df`).
    :return: List of :class:`FastaEntry` objects, one per row.
    :rtype: list[FastaEntry]
    """

    missing_cols = set(COLS) - set(df.columns)
    if missing_cols:
        message = (
            f"The following expected columns are missing from the dataframe: {missing_cols}. "
            "These will be filled with default values in the resulting FastaEntry objects."
        )
        logger.warning("df_to_entries: %s", message)
        warnings.warn(message, stacklevel=2)

    cols = [c for c in COLS if c in df.columns]
    entries = [FastaEntry(**row.to_dict()) for _, row in df[cols].iterrows()]
    logger.debug("df_to_entries: converted %d rows to FastaEntry objects", len(entries))
    return entries


def entries_to_fasta(
    entries: Iterable[FastaEntry],
    output_file: str | None = None,
    line_wrapping: bool = True,
    max_sequence_length: int | None = None,
) -> StringIO | None:
    """
    Serialize :class:`FastaEntry` objects to FASTA text or a file.

    :param entries: The :class:`FastaEntry` objects to serialize.
    :param output_file: Destination path. If None (default), the FASTA text is
        returned as a :class:`io.StringIO`.
    :param line_wrapping: Whether sequence wrapping is allowed. Default True.
    :param max_sequence_length: Line width to wrap sequences at; ``None`` or a
        non-positive value means no wrapping. See :meth:`FastaEntry.serialize`.
    :return: A :class:`io.StringIO` with the FASTA content, or None if
        ``output_file`` was provided.
    :rtype: StringIO | None
    :raises OSError: If ``output_file`` cannot be written.
    """

    fasta_string = StringIO()
    count = 0
    for entry in entries:
        fasta_string.write(entry.serialize(line_wrapping=line_wrapping, max_sequence_length=max_sequence_length))
        count += 1

    fasta_string.seek(0)

    if output_file is not None:
        with open(file=output_file, mode="w", encoding="UTF-8") as f:
            f.write(fasta_string.getvalue())
        logger.debug("entries_to_fasta: wrote %d entries to %r", count, output_file)
        return None

    logger.debug("entries_to_fasta: serialized %d entries to StringIO", count)
    return fasta_string


def to_fasta(
    data: pd.DataFrame | Iterable[FastaEntry],
    output_file: str | None = None,
    line_wrapping: bool = True,
    max_sequence_length: int | None = None,
) -> StringIO | None:
    """
    Convert a FASTA DataFrame or list of entries to FASTA text or a file.

    :param data: The FASTA DataFrame or an iterable of :class:`FastaEntry` objects.
    :param output_file: Destination path. If None, the FASTA text is returned as
        a :class:`io.StringIO`.
    :param line_wrapping: Whether sequence wrapping is allowed. Default True.
    :param max_sequence_length: Line width to wrap sequences at; ``None`` means
        no wrapping.
    :return: A :class:`io.StringIO` with the FASTA content, or None if
        ``output_file`` was provided.
    :rtype: StringIO | None
    """

    if isinstance(data, pd.DataFrame):
        return entries_to_fasta(df_to_entries(data), output_file, line_wrapping, max_sequence_length)

    return entries_to_fasta(data, output_file, line_wrapping, max_sequence_length)


def _extract_fasta_header_elements(entry_str: str) -> list[str]:
    """
    Split a header line into space-separated tokens (``>`` stripped).

    :param entry_str: The header line of a fasta entry.
    :return: List of elements extracted from the header line.
    :rtype: list[str]
    """

    line_elements = entry_str.rstrip().replace(">", "").split(" ")
    return line_elements


def _extract_initial_info(line_elements: list[str]) -> tuple[str | None, str, str | None]:
    """
    Extract ``(db, unique_identifier, entry_name)`` from the header elements.

    Expects the first element to be ``db|unique_identifier|entry_name``. If it
    does not split into exactly three parts, the whole first element is used as
    the ``unique_identifier`` and a :class:`UserWarning` is emitted.

    :param line_elements: Elements from :func:`_extract_fasta_header_elements`.
    :return: Tuple containing database, unique identifier, and entry name.
    :rtype: tuple[str | None, str, str | None]
    """

    first_element_parts = line_elements[0].split("|")

    if len(first_element_parts) == 3:
        db = first_element_parts[0]
        unique_identifier = first_element_parts[1]
        entry_name = first_element_parts[2]
        return db, unique_identifier, entry_name

    logger.warning("Non-standard FASTA header %r; using it as the unique_identifier", line_elements[0])
    warnings.warn(
        f"Invalid fasta header format: {line_elements[0]}, using only the first part as unique identifier.",
        stacklevel=2,
    )

    return None, line_elements[0], None


def _process_line_elements(line_elements: list[str]) -> dict[str, list[str]]:
    """
    Process the header tokens and group them by field key.

    The first element (``>db|UniqueIdentifier|EntryName``) is skipped. Tokens
    before the first ``KEY=`` are the protein name; ``KEY=value`` tokens start a
    new field. Keys outside the standard set are kept as-is so they can be
    preserved in :attr:`FastaEntry.additional_fields`.

    :param line_elements: Elements from :func:`_extract_fasta_header_elements`.
    :return: Dictionary mapping each field key to its list of value tokens.
    :rtype: dict[str, list[str]]
    """

    info: dict[str, list[str]] = {}
    current_state = FastaFields.PROTEIN_NAME.value

    for elem in line_elements[1:]:
        value = elem
        if "=" in elem:
            # Split on the first '=' so keys of any length round-trip.
            current_state, _, value = elem.partition("=")

        if current_state not in {f.value for f in FastaFields}:
            logger.debug("Unrecognized header field %r; captured in additional_fields", current_state)

        info.setdefault(current_state, []).append(value)

    return info


def _fasta_str_to_entry(fasta_str: str) -> FastaEntry:
    """
    Parse a single FASTA header line into a :class:`FastaEntry`.

    Standard keys (``PN``/``OS``/``OX``/``GN``/``PE``/``SV``) map to named
    attributes; any other key is preserved in
    :attr:`FastaEntry.additional_fields`.

    :param fasta_str: The header line of a fasta entry (leading ``>`` included).
    :return: FastaEntry object containing the extracted information.
    :rtype: FastaEntry
    :raises FastaFormatError: If the header has no usable identifier at all.
    """

    line_elements = _extract_fasta_header_elements(fasta_str)
    db, unique_identifier, entry_name = _extract_initial_info(line_elements)
    info = _process_line_elements(line_elements)

    if not db and not unique_identifier and not entry_name:
        raise FastaFormatError(fasta_str, reason="no db, unique_identifier, or entry_name found")

    def _join_list_values(data: dict[str, list[str]]) -> dict[str, str | None]:
        return {k: " ".join(v) if v else None for k, v in data.items()}

    joined_info = _join_list_values(info)
    standard_keys = {f.value for f in FastaFields}
    additional_fields = {k: v for k, v in joined_info.items() if k not in standard_keys and v}

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
        additional_fields=additional_fields,
    )
