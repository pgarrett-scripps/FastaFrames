"""
Core fastaframes functionality: the :class:`FastaEntry` record and the functions
that convert between FASTA text, ``FastaEntry`` objects, and pandas DataFrames.

Parsing and writing are done by `fastatacular <https://pypi.org/project/fastatacular/>`_;
this module maps its entries to the FastaFrames column names.

Errors
------
Malformed FASTA (an empty header, sequence data before the first header) raises
:class:`fastaframes.exceptions.FastaFormatError` (a subclass of ``ValueError``).
Unrecognized ``KEY=value`` header fields are kept in
:attr:`FastaEntry.additional_fields`. Writing an entry that cannot be written (no
identifier, a non-integer ``organism_identifier``, ...) raises
:class:`fastatacular.FastaWriteError` (also a ``ValueError``).

Logging
-------
Logs under the ``"fastaframes.fastaframes"`` logger: normal progress (entry
counts, output targets) at ``DEBUG`` and tolerated-but-odd input at ``WARNING``.
"""

import ast
import logging
import os
import warnings
from collections.abc import Generator, Iterable, Iterator, Mapping
from dataclasses import asdict, dataclass, field
from enum import Enum
from io import StringIO
from typing import Any

import pandas as pd
from fastatacular import FastaParseError, FastaReader, SequenceEntry, from_records, read_fasta, write_fasta

from .exceptions import FastaFormatError
from .util import as_source, get_lines

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

_STR_COLS = [c for c in COLS if c != "protein_sequence"]


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
    An identifier that is not exactly three ``|``-separated fields (``>myid``,
    ``>gi|229552|pdb|1MBA|A``) is kept whole in :attr:`unique_identifier`, with
    ``db`` and ``entry_name`` set to None. Unrecognized ``KEY=value`` header fields
    are kept in :attr:`additional_fields` and written back after the standard ones.

    :param db: Source database (``"sp"`` for Swiss-Prot, ``"tr"`` for TrEMBL).
    :param unique_identifier: Primary accession (e.g. ``"P12345"``).
    :param entry_name: Entry name (e.g. ``"CP2D7_HUMAN"``).
    :param protein_name: Recommended protein name (the text before the first ``KEY=``).
    :param organism_name: Scientific organism name (``OS=``).
    :param organism_identifier: NCBI taxonomy id (``OX=``).
    :param gene_name: Gene name (``GN=``).
    :param protein_existence: Protein-existence evidence level (``PE=``).
    :param sequence_version: Sequence version (``SV=``).
    :param protein_sequence: The amino-acid sequence.
    :param additional_fields: Non-standard ``KEY=value`` header fields. Defaults to
        an empty dict.
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
        """The FASTA header line for this entry (leading ``>`` included, no newline).

        The identifier is ``db|unique_identifier|entry_name`` when ``db`` or
        ``entry_name`` is set, else ``unique_identifier`` alone. Then come the
        protein name, ``OS=``, ``OX=``, ``GN=``, ``PE=``, ``SV=`` (UniProt order)
        and :attr:`additional_fields`; unset fields are left out.

        :return: The FASTA header line.
        :rtype: str
        :raises fastatacular.FastaWriteError: If the entry has no identifier.
        """
        return ">" + _to_sequence_entries([self])[0].raw_header

    def serialize(self, line_wrapping: bool = True, max_sequence_length: int | None = None) -> str:
        """Serialize this entry to a FASTA string: :attr:`header`, then the sequence.

        The sequence is wrapped only when ``line_wrapping`` is True and
        ``max_sequence_length`` is a positive int; otherwise it is written on one
        line. Every line ends with ``\\n``; an empty sequence writes the header alone.

        :param line_wrapping: Whether wrapping is allowed at all. Default True.
        :param max_sequence_length: Line width to wrap the sequence at. ``None``
            (default) or a non-positive value means "do not wrap".
        :return: The serialized FASTA record.
        :rtype: str
        """
        out = StringIO()
        write_fasta(_to_sequence_entries([self]), out, line_width=_line_width(line_wrapping, max_sequence_length))
        return out.getvalue()

    def to_dict(self) -> dict[str, object]:
        """Return this entry as a dict of all fields plus the computed ``protein_id``.

        :return: Dictionary representation of the entry (includes ``protein_id``).
        :rtype: dict[str, object]
        """
        d: dict[str, object] = asdict(self)
        d["protein_id"] = self.protein_id
        return d


# --- conversion to and from fastatacular -------------------------------------------------


def _is_missing(value: object) -> bool:
    """None, NaN and ``pd.NA`` are missing; containers are not."""
    if value is None:
        return True
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def _text(value: object) -> str | None:
    """A str for a cell value: missing -> None, ``9606.0`` -> ``"9606"``."""
    if _is_missing(value):
        return None
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def _fields(value: object) -> dict[str, str]:
    """``additional_fields`` from a dict, its ``str()`` (a CSV round trip), or ``KEY=value`` text."""
    if isinstance(value, Mapping):
        return {str(k): str(v) for k, v in value.items()}
    if _is_missing(value):
        return {}
    text = str(value).strip()
    if not text:
        return {}
    if text.startswith("{"):
        parsed = ast.literal_eval(text)
        if isinstance(parsed, dict):
            return {str(k): str(v) for k, v in parsed.items()}
    return dict(read_fasta(StringIO(f">x {text}\n"))[0].extra)


def _from_sequence_entry(entry: SequenceEntry) -> FastaEntry:
    parts = entry.identifier.split("|")
    db, uid, name = parts if len(parts) == 3 else (None, entry.identifier, None)
    extra = dict(entry.extra)
    protein_name = entry.pname
    if protein_name is None and "PN" in extra:
        # Headers written by FastaFrames 1.x: ">sp|P1|ABC PN=Alpha protein OS=...".
        protein_name = extra.pop("PN")
    return FastaEntry(
        db=db,
        unique_identifier=uid,
        entry_name=name,
        protein_name=protein_name,
        organism_name=entry.os_name,
        organism_identifier=_text(entry.ncbi_tax_id),
        gene_name=entry.gname,
        protein_existence=_text(entry.pe),
        sequence_version=_text(entry.sv),
        protein_sequence=entry.sequence,
        additional_fields=extra,
    )


def _to_record(entry: FastaEntry) -> dict[str, object]:
    db, uid, name = _text(entry.db), _text(entry.unique_identifier), _text(entry.entry_name)
    identifier = f"{db or ''}|{uid or ''}|{name or ''}" if (db or name) else uid
    return {
        "identifier": identifier,
        "pname": entry.protein_name,
        "os_name": entry.organism_name,
        "ncbi_tax_id": entry.organism_identifier,
        "gname": entry.gene_name,
        "pe": entry.protein_existence,
        "sv": entry.sequence_version,
        "extra": _fields(entry.additional_fields) or None,
        "sequence": entry.protein_sequence,
    }


def _to_sequence_entries(entries: Iterable[FastaEntry]) -> list[SequenceEntry]:
    return from_records(_to_record(e) for e in entries)


def _line_width(line_wrapping: bool, max_sequence_length: int | None) -> int:
    if line_wrapping and max_sequence_length is not None and max_sequence_length > 0:
        return int(max_sequence_length)
    return 0


def _format_error(err: FastaParseError) -> FastaFormatError:
    exc = FastaFormatError((err.context or "").rstrip("\r\n"), reason=str(err))
    exc.__cause__ = err
    return exc


def _without_bad_records(lines: Iterable[str]) -> Iterator[str]:
    """Drop headers with no identifier (and their sequence) and sequence lines before the first header."""
    keep = True
    seen_header = False
    for line in lines:
        if line.startswith(">"):
            seen_header = True
            keep = bool(line[1:].strip())
            if not keep:
                logger.warning("Skipping malformed FASTA header %r: no identifier", line)
        elif not seen_header and line.strip() and line[0] not in ";#":
            logger.warning("Skipping sequence line before the first header: %r", line)
            continue
        if keep:
            yield line + "\n"


def _fasta_str_to_entry(fasta_str: str) -> FastaEntry:
    """
    Parse a single FASTA header line into a :class:`FastaEntry`.

    :param fasta_str: The header line of a fasta entry (leading ``>`` included).
    :return: FastaEntry object with an empty sequence.
    :rtype: FastaEntry
    :raises FastaFormatError: If the header has no identifier.
    """
    try:
        entries = read_fasta(StringIO(fasta_str.rstrip("\r\n") + "\n"))
    except FastaParseError as err:
        raise _format_error(err) from err
    if len(entries) != 1:
        raise FastaFormatError(fasta_str, reason="expected exactly one header line")
    return _from_sequence_entry(entries[0])


# --- public API ------------------------------------------------------------------------


def fasta_to_entries(data: Any, skip_error: bool = False) -> Generator[FastaEntry, None, None]:
    """
    Parse FASTA content and lazily yield one :class:`FastaEntry` per record.

    :param data: A path (``str`` or :class:`os.PathLike`; plain, gzip, bzip2 or xz),
        raw FASTA text, a readable text or binary file-like object (read from the
        start), or an iterable of ``str``/``bytes`` lines. A ``str`` is FASTA text
        when it is blank, contains a newline or starts with ``>``; otherwise it is a path.
    :param skip_error: If True, a header with no identifier (and its sequence) and
        sequence lines before the first header are skipped with a warning; if False
        (default) they raise.
    :return: A generator that yields :class:`FastaEntry` objects.
    :rtype: Generator[FastaEntry, None, None]
    :raises FastaFormatError: On malformed FASTA when ``skip_error`` is False.
    :raises FileNotFoundError: If ``data`` is a path that does not exist.
    :raises InvalidInputError: If ``data`` is not a supported source.
    """
    source = StringIO("".join(_without_bad_records(get_lines(data)))) if skip_error else as_source(data)
    yielded = 0
    try:
        with FastaReader(source) as reader:
            for entry in reader:
                yield _from_sequence_entry(entry)
                yielded += 1
    except FastaParseError as err:
        raise _format_error(err) from err
    logger.debug("fasta_to_entries: yielded %d entries", yielded)


def entries_to_df(entries: Iterable[FastaEntry]) -> pd.DataFrame:
    """
    Build a DataFrame from :class:`FastaEntry` objects (one row per entry).

    Columns are the :data:`COLS` fields plus ``additional_fields`` and the
    computed ``protein_id`` (see :meth:`FastaEntry.to_dict`).

    :param entries: The :class:`FastaEntry` objects to tabulate.
    :return: A pandas DataFrame representing the FastaEntry objects.
    :rtype: pd.DataFrame
    """
    return pd.DataFrame([entry.to_dict() for entry in entries])


def to_df(data: Any, skip_error: bool = False) -> pd.DataFrame:
    """
    Convert a FASTA source (or list of entries) into a pandas DataFrame.

    :param data: Any source :func:`fasta_to_entries` accepts, or a list of
        :class:`FastaEntry` objects.
    :param skip_error: If True, malformed FASTA records are skipped instead of
        raising. Ignored when ``data`` is already a list of entries.
    :return: A pandas DataFrame representing the FASTA content or entries.
    :rtype: pd.DataFrame
    :raises FastaFormatError: On malformed FASTA when ``skip_error`` is False.
    :raises FileNotFoundError: If ``data`` is a path that does not exist.
    """
    if isinstance(data, list):
        return entries_to_df(data)
    return entries_to_df(fasta_to_entries(data, skip_error))


def df_to_entries(df: pd.DataFrame) -> list[FastaEntry]:
    """
    Convert a FASTA DataFrame into a list of :class:`FastaEntry` objects.

    The :data:`COLS` columns and ``additional_fields`` are used; other columns
    (such as ``protein_id``) are ignored. Missing values (None, NaN, ``pd.NA``)
    become None (an empty sequence for ``protein_sequence``), and integral
    floats such as ``9606.0`` become ``"9606"``. ``additional_fields`` may be a
    dict, its ``str()`` (as read back from CSV) or ``KEY=value`` text. Missing
    expected columns are left at their defaults with a :class:`UserWarning`.

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

    entries = []
    for row in df.to_dict("records"):
        kwargs: dict[str, Any] = {c: _text(row[c]) for c in _STR_COLS if c in row}
        if "unique_identifier" in kwargs and kwargs["unique_identifier"] is None:
            kwargs["unique_identifier"] = ""
        if "protein_sequence" in row:
            kwargs["protein_sequence"] = _text(row["protein_sequence"]) or ""
        if "additional_fields" in row:
            kwargs["additional_fields"] = _fields(row["additional_fields"])
        entries.append(FastaEntry(**kwargs))
    logger.debug("df_to_entries: converted %d rows to FastaEntry objects", len(entries))
    return entries


def entries_to_fasta(
    entries: Iterable[FastaEntry],
    output_file: str | os.PathLike[str] | None = None,
    line_wrapping: bool = True,
    max_sequence_length: int | None = None,
) -> StringIO | None:
    """
    Serialize :class:`FastaEntry` objects to FASTA text or a file.

    Every entry is checked before anything is written.

    :param entries: The :class:`FastaEntry` objects to serialize.
    :param output_file: Destination path (``.gz``, ``.bz2`` and ``.xz`` are
        compressed). If None (default), the FASTA text is returned as a
        :class:`io.StringIO`.
    :param line_wrapping: Whether sequence wrapping is allowed. Default True.
    :param max_sequence_length: Line width to wrap sequences at; ``None`` or a
        non-positive value means no wrapping. See :meth:`FastaEntry.serialize`.
    :return: A :class:`io.StringIO` with the FASTA content, or None if
        ``output_file`` was provided.
    :rtype: StringIO | None
    :raises fastatacular.FastaWriteError: If an entry cannot be written.
    :raises OSError: If ``output_file`` cannot be written.
    """
    seq_entries = _to_sequence_entries(entries)
    width = _line_width(line_wrapping, max_sequence_length)

    if output_file is not None:
        write_fasta(seq_entries, os.fspath(output_file), line_width=width)
        logger.debug("entries_to_fasta: wrote %d entries to %r", len(seq_entries), output_file)
        return None

    fasta_string = StringIO()
    write_fasta(seq_entries, fasta_string, line_width=width)
    fasta_string.seek(0)
    logger.debug("entries_to_fasta: serialized %d entries to StringIO", len(seq_entries))
    return fasta_string


def to_fasta(
    data: pd.DataFrame | Iterable[FastaEntry],
    output_file: str | os.PathLike[str] | None = None,
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
