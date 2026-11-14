"""A short tour of the fastaframes API.

Run it from anywhere: the input is resolved relative to this file, and the output is
written to a temporary directory so that running the example never leaves files in the
repository.
"""

import tempfile
from io import StringIO, TextIOWrapper
from pathlib import Path

from fastaframes import to_df, to_fasta

# Resolve the input relative to this file, so the example runs from any working directory.
FASTA_FILE = Path(__file__).resolve().parent.parent / "tests" / "data" / "test.fasta"
OUTPUT_DIR = Path(tempfile.mkdtemp(prefix="fastaframes-tutorial-"))

# fastaframes accepts multiple types of inputs, in this example we will use a file path
fasta_df = to_df(data=FASTA_FILE)

# It also supports IO like inputs
with open(FASTA_FILE) as f:
    assert isinstance(f, TextIOWrapper)
    fasta_df2 = to_df(data=f)
assert fasta_df.equals(fasta_df2)

# Here are the possible columns:
cols = [
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
    "additional_fields",
    "raw_header",
    "protein_id",
]
assert list(fasta_df.columns) == cols

ids = fasta_df["unique_identifier"].tolist()
assert ids == ["A0A087X1C5", "A0A0B4J2F2", "A0A0C5B5G6", "A0A0K2S4Q6"]

# Directly write to a fasta file
output_file = OUTPUT_DIR / "test.fasta"
to_fasta(fasta_df, output_file)

# or create an IO object
fasta_io = to_fasta(fasta_df)
assert isinstance(fasta_io, StringIO)

print(f"{len(fasta_df)} entries, columns {len(cols)}, written to {output_file}")
