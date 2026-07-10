"""Tests for the 1.4.0 additions: additional_fields, line wrapping, exceptions, logging."""

import io
import logging

import pytest

from fastaframes import (
    FastaEntry,
    FastaFormatError,
    FastaFramesError,
    InvalidInputError,
    fasta_to_entries,
    to_fasta,
)
from fastaframes.fastaframes import _fasta_str_to_entry
from fastaframes.util import get_lines


def test_additional_fields_round_trip():
    """Non-standard KEY=value header fields survive serialize -> parse."""
    entry = FastaEntry(
        db="sp",
        unique_identifier="A",
        entry_name="B",
        protein_sequence="SEQ",
        additional_fields={"XT": "kinase"},
    )
    serialized = entry.serialize()
    assert "XT=kinase" in serialized

    reparsed = next(fasta_to_entries(io.StringIO(serialized)))
    assert reparsed.additional_fields == {"XT": "kinase"}


def test_unknown_field_is_stored_not_raised():
    """A non-standard field no longer raises; it lands in additional_fields."""
    entry = _fasta_str_to_entry(">sp|A|B Name ZZ=weird")
    assert entry.additional_fields == {"ZZ": "weird"}


def test_line_wrapping_wraps_at_width():
    entry = FastaEntry(db="sp", unique_identifier="A", entry_name="B", protein_sequence="ABCDEFGHIJ")
    out = entry.serialize(max_sequence_length=4)
    body = out.splitlines()[1:]
    assert body == ["ABCD", "EFGH", "IJ"]


def test_line_wrapping_disabled_by_default_and_by_zero():
    entry = FastaEntry(db="sp", unique_identifier="A", entry_name="B", protein_sequence="ABCDEFGHIJ")
    assert entry.serialize().splitlines()[1] == "ABCDEFGHIJ"
    assert entry.serialize(max_sequence_length=0).splitlines()[1] == "ABCDEFGHIJ"


def test_to_fasta_wrapping_multiple_entries():
    e1 = FastaEntry(db="sp", unique_identifier="A", entry_name="B", protein_sequence="AAAAAA")
    e2 = FastaEntry(db="sp", unique_identifier="C", entry_name="D", protein_sequence="CCCCCC")
    out = to_fasta([e1, e2], max_sequence_length=3).getvalue()
    assert out == ">sp|A|B\nAAA\nAAA\n>sp|C|D\nCCC\nCCC\n"


def test_header_does_not_emit_literal_none():
    entry = FastaEntry(db=None, unique_identifier="Q8I6R7", entry_name=None)
    assert entry.header == ">|Q8I6R7|"
    assert "None" not in entry.serialize()


def test_exceptions_are_importable_and_subclassed():
    assert issubclass(FastaFormatError, FastaFramesError)
    assert issubclass(FastaFormatError, ValueError)
    assert issubclass(InvalidInputError, FastaFramesError)
    assert issubclass(InvalidInputError, TypeError)


def test_fasta_format_error_carries_attributes():
    with pytest.raises(FastaFormatError) as exc_info:
        _fasta_str_to_entry(">")
    assert exc_info.value.reason


def test_invalid_input_error_on_unsupported_type():
    with pytest.raises(InvalidInputError) as exc_info:
        list(get_lines(42))
    assert exc_info.value.received_type is int


def test_get_lines_decodes_binary_stream():
    assert list(get_lines(io.BytesIO(b">sp|A|B\nSEQ\n"))) == [">sp|A|B", "SEQ"]


def test_skip_error_swallows_format_error():
    df_entries = list(fasta_to_entries(io.StringIO(">\nSEQ\n>sp|A|B\nSEQ\n"), skip_error=True))
    assert [e.unique_identifier for e in df_entries] == ["A"]


def test_logging_emitted_when_enabled(caplog):
    with caplog.at_level(logging.DEBUG, logger="fastaframes"):
        list(fasta_to_entries(io.StringIO(">sp|A|B\nSEQ\n")))
    assert any("yielded" in rec.getMessage() for rec in caplog.records)
