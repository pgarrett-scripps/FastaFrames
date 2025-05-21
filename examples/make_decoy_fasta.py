import pandas as pd

from fastaframes import to_df, to_fasta

FASTA_FILE = "../tests/data/test.fasta"

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
to_fasta(df, "decoy.fasta")
