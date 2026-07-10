![CI](https://github.com/pgarrett-scripps/FastaFrames/actions/workflows/ci.yml/badge.svg)
[![PyPI](https://img.shields.io/pypi/v/fastaframes)](https://pypi.org/project/fastaframes/)
[![Python](https://img.shields.io/pypi/pyversions/fastaframes)](https://pypi.org/project/fastaframes/)

# FastaFrames

Convert between UniProt FASTA files and pandas DataFrames.

## Installation

```sh
pip install fastaframes
```

## Quick Start

### Read a FASTA file into a DataFrame

```python
from fastaframes import to_df

df = to_df("proteins.fasta")
print(df.head())
```

### Write a DataFrame back to FASTA

```python
from fastaframes import to_fasta

to_fasta(df, output_file="output.fasta")
```

### Work with individual entries

```python
from fastaframes import fasta_to_entries, entries_to_fasta

for entry in fasta_to_entries("proteins.fasta"):
    print(entry.unique_identifier, entry.protein_name)

# Filter and write back
entries = [e for e in fasta_to_entries("proteins.fasta") if e.organism_name == "Homo sapiens"]
entries_to_fasta(entries, output_file="human_only.fasta")
```

### Multiple input formats

```python
from io import StringIO
from fastaframes import to_df

# From a file path
df = to_df("proteins.fasta")

# From a string
df = to_df(">sp|P12345|EXAMPLE_HUMAN Example protein OS=Homo sapiens OX=9606\nMSEQUENCE\n")

# From a file object
with open("proteins.fasta") as f:
    df = to_df(f)

# From a StringIO
df = to_df(StringIO(">sp|P12345|EXAMPLE_HUMAN\nMSEQUENCE\n"))
```

### Skip malformed entries

```python
from fastaframes import to_df

df = to_df("messy_data.fasta", skip_error=True)
```

### Wrap sequence lines

```python
from fastaframes import to_fasta

# Wrap sequences at 60 characters (default is a single line per record)
to_fasta(df, output_file="wrapped.fasta", max_sequence_length=60)
```

### Preserve non-standard header fields

Any `KEY=value` field beyond the standard `PN`/`OS`/`OX`/`GN`/`PE`/`SV` is kept
in `FastaEntry.additional_fields` and re-emitted on serialization.

### Error handling

All package errors derive from `FastaFramesError`. Malformed FASTA raises
`FastaFormatError` (also a `ValueError`); unsupported inputs raise
`InvalidInputError` (also a `TypeError`).

```python
from fastaframes import to_df, FastaFormatError

try:
    df = to_df("example.fasta")
except FastaFormatError as err:
    print(err.header, err.reason)
```

### Logging

FastaFrames logs under the `fastaframes` logger and stays silent until you
configure logging:

```python
import logging
logging.getLogger("fastaframes").setLevel(logging.DEBUG)
```

## DataFrame Columns

Given this FASTA entry:

```
>sp|A0A087X1C5|CP2D7_HUMAN Putative cytochrome P450 2D7 OS=Homo sapiens OX=9606 GN=CYP2D7 PE=5 SV=1
MGLEALVPLAMIVAIFLLLVDLMHRHQRWAARYPPGPLPLPGLGNLLHVDFQNTPYCFDQ
```

`to_df` produces:

| db | unique_identifier | entry_name  | protein_name                 | organism_name | organism_identifier | gene_name | protein_existence | sequence_version | protein_sequence                                                 |
|----|-------------------|-------------|------------------------------|---------------|---------------------|-----------|-------------------|------------------|------------------------------------------------------------------|
| sp | A0A087X1C5        | CP2D7_HUMAN | Putative cytochrome P450 2D7 | Homo sapiens  | 9606                | CYP2D7    | 5                 | 1                | MGLEALVPLAMIVAIFLLLVDLMHRHQRWAARYPPGPLPLPGLGNLLHVDFQNTPYCFDQ |

Column descriptions (following the [UniProt FASTA header format](https://www.uniprot.org/help/fasta-headers)):

| Column | Description |
|--------|-------------|
| `db` | Database source: `sp` (Swiss-Prot) or `tr` (TrEMBL) |
| `unique_identifier` | Primary UniProtKB accession number |
| `entry_name` | UniProtKB entry name |
| `protein_name` | Recommended protein name (RecName or first SubName) |
| `organism_name` | Scientific name of the source organism |
| `organism_identifier` | NCBI taxonomy identifier |
| `gene_name` | First gene name (if available) |
| `protein_existence` | Numerical evidence code for protein existence |
| `sequence_version` | Sequence version number |
| `protein_sequence` | Amino acid sequence |

## Development

```sh
pip install -e ".[dev]"
```

Common commands via [just](https://github.com/casey/just):

```sh
just check      # Run all checks (lint, typecheck, test)
just lint       # Lint with ruff
just fmt        # Format with ruff
just typecheck  # Type check with ty
just test       # Run tests
just test -v    # Run tests verbosely
```
