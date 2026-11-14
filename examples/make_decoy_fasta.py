"""Build a target-decoy FASTA database from a target one.

Run it from anywhere: the input is resolved relative to this file, and the output is
written to a temporary directory so that running the example never leaves files in the
repository.
"""

import tempfile
from pathlib import Path

import pandas as pd

from fastaframes import to_df, to_fasta

# Resolve the input relative to this file, so the example runs from any working directory.
FASTA_FILE = Path(__file__).resolve().parent.parent / "tests" / "data" / "test.fasta"
OUTPUT_DIR = Path(tempfile.mkdtemp(prefix="fastaframes-decoy-"))

# fastaframes accepts multiple types of inputs, in this example we will use a file path
fasta_df = to_df(data=FASTA_FILE)

# Copy the dataframe to create a decoy database
decoy_df = fasta_df.copy(deep=True)

# Append "decoy_" to the start of db
decoy_df["db"] = "decoy_" + decoy_df["db"]

# Reverse the protein sequences
decoy_df["protein_sequence"] = decoy_df["protein_sequence"].str[::-1]

# Concatenate the two dataframes
df = pd.concat([fasta_df, decoy_df])

# write df to FASTA file
output_file = OUTPUT_DIR / "decoy.fasta"
to_fasta(df, output_file)

print(f"{len(fasta_df)} target and {len(decoy_df)} decoy entries written to {output_file}")
