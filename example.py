from fastaframes import to_df, to_fasta

# Read IO input
with open('example.fasta', 'r') as fasta_io:
    fasta_df = to_df(fasta_data=fasta_io)

# Read File input
fasta_df = to_df(fasta_data='example.fasta')

fasta_df.to_csv('output.csv')  # output fasta in csv format
print(fasta_df.head())

# Write StringIO to file
fasta_io = to_fasta(fasta_data=fasta_df)  # outputs StringIO if file=None
with open('output.fasta', 'w') as output_file:
    output_file.write(fasta_io.getvalue())

# Write directly to file
to_fasta(fasta_data=fasta_df, file='output.fasta')
