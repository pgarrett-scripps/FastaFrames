![CI](https://github.com/pgarrett-scripps/FastaFrames/actions/workflows/ci.yml/badge.svg)
[![PyPI](https://img.shields.io/pypi/v/fastaframes)](https://pypi.org/project/fastaframes/)
[![Python](https://img.shields.io/pypi/pyversions/fastaframes)](https://pypi.org/project/fastaframes/)

# FastaFrames

FastaFrames reads protein FASTA files into pandas DataFrames and writes
DataFrames back to FASTA. Each record becomes one row, and each field of a
[UniProt FASTA header](https://www.uniprot.org/help/fasta-headers) (accession,
entry name, protein name, organism, gene name, protein existence, sequence
version) becomes one column. Parsing and writing are done by
[fastatacular](https://pypi.org/project/fastatacular/); FastaFrames maps its
entries to DataFrame columns.

FastaFrames requires Python 3.12 or later.

## Installation

```sh
pip install fastaframes
```

## Quick start

The examples below use this file, saved as `proteins.fasta`:

```
>sp|A0A087X1C5|CP2D7_HUMAN Putative cytochrome P450 2D7 OS=Homo sapiens OX=9606 GN=CYP2D7 PE=5 SV=1
MGLEALVPLAMIVAIFLLLVDLMHRHQRWAARYPPGPLPLPGLGNLLHVDFQNTPYCFDQ
>sp|A0A0C5B5G6|MOTSC_HUMAN Mitochondrial-derived peptide MOTS-c OS=Homo sapiens OX=9606 GN=MT-RNR1 PE=1 SV=1
MRWQEMGYIFYPRKLR
```

### Read a FASTA file into a DataFrame

```python
from fastaframes import to_df

df = to_df("proteins.fasta")
print(df[["unique_identifier", "gene_name", "protein_existence"]].to_string())
```

```
  unique_identifier gene_name protein_existence
0        A0A087X1C5    CYP2D7                 5
1        A0A0C5B5G6   MT-RNR1                 1
```

`to_df` also accepts a `pathlib.Path`, a gzip, bzip2 or xz compressed file, an
open file object, or FASTA text (a string that starts with `>` or contains a
newline). Any other string is treated as a path, and a missing path raises
`FileNotFoundError`.

### Write a DataFrame to FASTA

```python
from fastaframes import to_fasta

short = df[df["protein_sequence"].str.len() < 30]
to_fasta(short, output_file="short.fasta")
print(open("short.fasta").read(), end="")
```

```
>sp|A0A0C5B5G6|MOTSC_HUMAN Mitochondrial-derived peptide MOTS-c OS=Homo sapiens OX=9606 GN=MT-RNR1 PE=1 SV=1
MRWQEMGYIFYPRKLR
```

Without `output_file`, `to_fasta` returns the FASTA text as an `io.StringIO`. An `output_file`
ending in `.gz`, `.bz2` or `.xz` is written compressed. Protein sequences are
written on one line unless `max_sequence_length` is set (for example
`max_sequence_length=60`).

### Round trip

UniProt FASTA files are written back byte for byte:

```python
from fastaframes import to_df, to_fasta

text = open("proteins.fasta").read()
assert to_fasta(to_df("proteins.fasta")).getvalue() == text
```

### Edit a header field

An edited column is changed in place in the original header; the other fields
keep their order and spacing:

```python
df = to_df("proteins.fasta")
df.loc[0, "gene_name"] = "CYP2D7P"
print(to_fasta(df).getvalue().splitlines()[0])
```

```
>sp|A0A087X1C5|CP2D7_HUMAN Putative cytochrome P450 2D7 OS=Homo sapiens OX=9606 GN=CYP2D7P PE=5 SV=1
```

### Work with entries

`fasta_to_entries` yields one `FastaEntry` dataclass per record, and
`entries_to_fasta` writes them back:

```python
from fastaframes import entries_to_fasta, fasta_to_entries

for entry in fasta_to_entries("proteins.fasta"):
    print(entry.unique_identifier, entry.protein_name)

entries = [e for e in fasta_to_entries("proteins.fasta") if e.protein_existence == "1"]
entries_to_fasta(entries, output_file="evidence_at_protein_level.fasta")
```

```
A0A087X1C5 Putative cytochrome P450 2D7
A0A0C5B5G6 Mitochondrial-derived peptide MOTS-c
```

### Errors

All package errors derive from `FastaFramesError`. Malformed FASTA raises
`FastaFormatError` (also a `ValueError`); `to_df(..., skip_error=True)` skips
malformed records with a warning instead. Unsupported inputs raise
`InvalidInputError` (also a `TypeError`). Writing an entry with no accession,
or a non-integer `organism_identifier`, `protein_existence` or
`sequence_version`, raises `fastatacular.FastaWriteError` (a `ValueError`).

FastaFrames logs under the `fastaframes` logger and is silent until logging is
configured.

## DataFrame columns

Column names follow the
[UniProt FASTA header format](https://www.uniprot.org/help/fasta-headers):

| Column | Content |
|--------|---------|
| `db` | Database: `sp` (UniProtKB/Swiss-Prot) or `tr` (UniProtKB/TrEMBL) |
| `unique_identifier` | UniProtKB accession, or the whole identifier for a non-UniProt header |
| `entry_name` | UniProtKB entry name |
| `protein_name` | Protein name |
| `organism_name` | Scientific name of the organism (`OS`) |
| `organism_identifier` | NCBI taxonomy identifier (`OX`) |
| `gene_name` | Gene name (`GN`) |
| `protein_existence` | Protein existence level, 1 to 5 (`PE`) |
| `sequence_version` | Sequence version (`SV`) |
| `protein_sequence` | Protein sequence |
| `additional_fields` | Other `KEY=value` header fields, as a dict |
| `raw_header` | The header as read, without `>` |
| `protein_id` | `db\|unique_identifier\|entry_name`, computed |

Header field values are strings (`organism_identifier` is `"9606"`), and
missing fields are `None`. `additional_fields` is written back after the
standard fields, also after a DataFrame or CSV round trip.

`raw_header` is used when writing: a row whose columns were not edited is
written as `raw_header`; an edited column is changed in place in it, and a
column set to `None` or `""` removes that field. To build headers from the
columns alone, drop `raw_header` or set it to `""`. DataFrames without a
`raw_header` column (for example from FastaFrames 1.x) are written from the
columns.

## Changes in 2.0

FastaFrames 2.0 parses and writes through fastatacular. The functions and
column names are the same as in 1.x, with one added column, `raw_header`. The
main differences:

1. Python 3.12 or later is required.
2. The protein name is no longer written as a `PN=` field, and a non-UniProt
   header is written as it was read (`>myid some description`). Files written
   by 1.x still read correctly.
3. A string that is not FASTA text is a path; a missing path raises
   `FileNotFoundError` instead of returning an empty DataFrame.
4. Sequence lines before the first header raise `FastaFormatError`.

The full list is in [CHANGELOG.md](CHANGELOG.md).

## Development

```sh
pip install -e ".[dev]"
just check   # ruff, ty and pytest
```

## License

MIT. See [LICENSE](LICENSE).

## Citation

If you use FastaFrames in published work, please cite it using the metadata in
[CITATION.cff](CITATION.cff).
