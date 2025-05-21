from io import TextIOWrapper, StringIO

from fastaframes import to_df
from fastaframes import to_fasta

FASTA_FILE = "../tests/data/test.fasta"

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
]

ids = fasta_df["unique_identifier"].tolist()
assert ids == ["A0A087X1C5", "A0A0B4J2F2", "A0A0C5B5G6", "A0A0K2S4Q6"]

# Directly write to a fasta file
to_fasta(fasta_df, "test.fasta")

# or create an IO object
fasta_io = to_fasta(fasta_df)
assert isinstance(fasta_io, StringIO)
