"""Regression tests for bugs in 1.4.0, fixed in 2.0.0 by parsing and writing through fastatacular."""

import gzip
import io
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from fastaframes import FastaEntry, FastaFormatError, fasta_to_entries, to_df, to_fasta

SP = ">sp|P12345|ABC_HUMAN Alpha protein OS=Homo sapiens OX=9606 GN=ABC PE=1 SV=2\nMKV\n"
SP2 = ">sp|P67890|DEF_HUMAN Beta protein OS=Homo sapiens OX=9606 PE=1 SV=1\nMKL\n"
UNIPROT_FILE = Path(__file__).parent / "data" / "test.fasta"


def _round_trip(text: str) -> str:
    return to_fasta(to_df(text)).getvalue()


def test_uniprot_headers_round_trip_byte_identical():
    """Bug 1: the protein name was written as an invented PN= field."""
    assert _round_trip(SP + SP2) == SP + SP2
    text = UNIPROT_FILE.read_text().rstrip("\n") + "\n"  # the file lacks a final newline
    assert to_fasta(to_df(UNIPROT_FILE), max_sequence_length=60).getvalue() == text


def test_pipeless_header_kept():
    """Bug 2: '>myid some desc' was rewritten as '>|myid| PN=some desc'."""
    assert _round_trip(">myid some desc\nMK\n") == ">myid some desc\nMK\n"


def test_missing_values_not_written_as_nan():
    """Bug 3: a column with gaps wrote GN=nan; a missing db/entry_name wrote '>nan|myid|nan'."""
    assert _round_trip(SP + SP2) == SP + SP2  # SP2 has no GN, so the gene_name column has a NaN
    df = pd.DataFrame({"db": [np.nan], "unique_identifier": ["myid"], "entry_name": [None], "protein_sequence": ["MK"]})
    with pytest.warns(UserWarning, match="missing"):
        assert to_fasta(df).getvalue() == ">myid\nMK\n"


def test_csv_round_trip():
    """Bug 3: after a CSV round trip OX=9606 came back as OX=9606.0 and gaps as nan."""
    buf = io.StringIO()
    to_df(SP + SP2 + ">sp|P1|X_HUMAN Gamma XX=1\nMK\n").to_csv(buf, index=False)
    buf.seek(0)
    assert to_fasta(pd.read_csv(buf)).getvalue() == SP + SP2 + ">sp|P1|X_HUMAN Gamma XX=1\nMK\n"


def test_additional_fields_survive_dataframe():
    """Bug 4: additional_fields were dropped by df_to_entries."""
    text = ">sp|P1|A_H Prot OS=Hs XX=1 Tag=a b\nMK\n"
    df = to_df(text)
    assert df.additional_fields[0] == {"XX": "1", "Tag": "a b"}
    assert to_fasta(df).getvalue() == text


def test_gt_inside_header_kept():
    """Bug 5: every '>' inside a header was deleted."""
    text = ">sp|P1|A_H Na(+)/H(+)->exchanger OS=Homo\nMK\n"
    assert to_df(text).protein_name[0] == "Na(+)/H(+)->exchanger"
    assert _round_trip(text) == text


def test_tab_separates_identifier():
    """Bug 6: a tab after the identifier was not a separator."""
    row = to_df(">sp|P1|A_H\tProt OS=Hs\nMK\n").iloc[0]
    assert (row.db, row.unique_identifier, row.entry_name, row.protein_name) == ("sp", "P1", "A_H", "Prot")


@pytest.mark.parametrize("as_path", [str, Path])
def test_gzip_and_pathlib_paths(tmp_path, as_path):
    """Bug 7: a .gz path crashed and pathlib.Path was rejected."""
    gz = tmp_path / "x.fasta.gz"
    with gzip.open(gz, "wt") as f:
        f.write(SP)
    plain = tmp_path / "x.fasta"
    plain.write_text(SP)
    assert to_df(as_path(gz)).unique_identifier.tolist() == ["P12345"]
    assert to_df(as_path(plain)).unique_identifier.tolist() == ["P12345"]


def test_missing_path_raises():
    """Bug 8: a missing path was read as FASTA text and gave an empty DataFrame."""
    with pytest.raises(FileNotFoundError):
        to_df("missing.fasta")
    with pytest.raises(FileNotFoundError):
        to_df(Path("missing.fasta"))


def test_nan_sequence_with_wrapping():
    """Bug 9: a NaN sequence with max_sequence_length raised TypeError."""
    df = pd.DataFrame({"unique_identifier": ["X", "Y"], "protein_sequence": [np.nan, "MKVL"]})
    with pytest.warns(UserWarning):
        assert to_fasta(df, max_sequence_length=2).getvalue() == ">X\n>Y\nMK\nVL\n"


def test_sequence_before_first_header():
    """Bug 9: sequence lines before the first header were silently dropped."""
    with pytest.raises(FastaFormatError):
        to_df("MKV\n>sp|P1|A x\nMK\n")
    assert to_df("MKV\n>sp|P1|A x\nMK\n", skip_error=True).unique_identifier.tolist() == ["P1"]


def test_reads_fastaframes_1x_output():
    """Files written by 1.x (PN= before the name) read back with the name in protein_name."""
    entry = next(fasta_to_entries(">sp|P1|A_H PN=Alpha protein OS=Hs\nMK\n"))
    assert entry.protein_name == "Alpha protein"
    assert entry.additional_fields == {}


def test_empty_sequence_round_trip():
    entries = [FastaEntry(db="sp", unique_identifier="A", entry_name="B"), FastaEntry(unique_identifier="C")]
    text = to_fasta(entries).getvalue()
    assert text == ">sp|A|B\n>C\n"
    assert list(fasta_to_entries(text)) == [
        FastaEntry(db="sp", unique_identifier="A", entry_name="B"),
        FastaEntry(db=None, unique_identifier="C", entry_name=None),
    ]


def test_raw_header_keeps_unedited_header_and_edits_in_place():
    """2.0 review: headers were rebuilt from columns, losing spacing and field order."""
    text = ">sp|P1|A_H  Prot  GN=G OS=Hs XX=1\nMK\n>sp|P2|B_H PN=Old name OS=Hs\nMK\n"
    df = to_df(text)
    assert df.raw_header[0] == "sp|P1|A_H  Prot  GN=G OS=Hs XX=1"
    assert to_fasta(df).getvalue() == text
    df.loc[0, "gene_name"] = "NEW"
    assert to_fasta(df).getvalue().startswith(">sp|P1|A_H  Prot  GN=NEW OS=Hs XX=1\n")


def test_malformed_additional_fields_cell_raises_format_error():
    df = to_df(SP)
    df.loc[0, "additional_fields"] = "{bad"
    with pytest.raises(FastaFormatError, match="row 0"):
        to_fasta(df)


def test_empty_gene_name_not_written():
    assert FastaEntry(unique_identifier="X", gene_name="").header == ">X"
    df = pd.DataFrame([FastaEntry(unique_identifier="X", gene_name="", protein_sequence="MK").to_dict()])
    assert to_fasta(df).getvalue() == ">X\nMK\n"


def test_leading_whitespace_before_header():
    """1.4.0 accepted indented headers; an indented one mid-file was merged into the sequence."""
    text = "  >sp|P1|A x\nMK\n >sp|P2|B y\n  MK\n"
    assert to_df(text).unique_identifier.tolist() == ["P1", "P2"]
    assert to_df(text, skip_error=True).protein_sequence.tolist() == ["MK", "MK"]


def test_repeated_key_keeps_first_and_round_trips():
    text = ">sp|P1|A X OS=H OS=M\nMK\n"
    assert to_df(text).organism_name[0] == "H"
    assert _round_trip(text) == text
