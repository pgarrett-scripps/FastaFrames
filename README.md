![example workflow](https://github.com/pgarrett-scripps/FastaFrames/actions/workflows/python-package.yml/badge.svg)
![example workflow](https://github.com/pgarrett-scripps/FastaFrames/actions/workflows/pylint.yml/badge.svg)

# FastaFrames
This Python module provides a set of functions to work with FASTA files. 
It allows you to read FASTA files, convert them to pandas dataframes, manipulate data, 
and write data back to FASTA files. It also supports converting FASTA files to a list of FastaEntry dataclass objects.

## Features
- Read FASTA files into pandas DataFrames
- Write FASTA files from pandas DataFrames

## Usage

To install fastaframes use pip:

```sh
pip install fastaframes
```

### Reading FASTA files
To read a FASTA file and convert it to a pandas DataFrame:

```python
from fastaframes import to_df

# IO input
with open('example.fasta', 'r') as fasta_io:
    fasta_df = to_df(fasta_data=fasta_io)

# or

# File input
fasta_df = to_df(fasta_data='example.fasta')
    
print(fasta_df.head())
```

### Writing FASTA files
To write a pandas DataFrame to a FASTA file:

```python
from fastaframes import to_fasta

# Write StringIO to file
fasta_io = to_fasta(fasta_data=fasta_df) # outputs StringIO if file=None
with open('output.fasta', 'w') as output_file:
    output_file.write(fasta_io.getvalue())

# or
    
# Write directly to file
to_fasta(fasta_data=fasta_df, file='output.fasta')
```

## Example DataFrame:

|   | db | unique_identifier | entry_name   | protein_name                                         | organism_name | organism_identifier | gene_name | protein_existence | sequence_version | protein_sequence                                       |
|---|----|------------------|--------------|------------------------------------------------------|---------------|---------------------|-----------|-------------------|------------------|--------------------------------------------------------|
| 0 | sp | A0A087X1C5       | CP2D7_HUMAN  | Putative cytochrome P450 2D7                         | Homo sapiens  | 9606.0              | CYP2D7    | 5.0               | 1.0              | MGLEALVPLAMIVAIFLLLVDLMHRHQRWAARYPPGPLPLPGLGNLLHVDFQNTPYCFDQ |
| 1 | sp | A0A0B4J2F2       | SIK1B_HUMAN  | Putative serine/threonine-protein kinase SIK1B        | Homo sapiens  | 9606.0              | SIK1B     | 5.0               | 1.0              | MVIMSEFSADPAGQGQGQQKPLRVGFYDIERTLGKGNFAVVKLARHRVTKTQVAIKIIDKLVQ |
| 2 | sp | A0A0C5B5G6       | MOTSC_HUMAN  | Mitochondrial-derived peptide MOTS-c                 | Homo sapiens  | 9606.0              | MT-RNR1   | 1.0               | 1.0              | MRWQEMGYIFYPRKLR                                      |
| 3 | sp | A0A0K2S4Q6       | CD3CH_HUMAN  | Protein CD300H                                       | Homo sapiens  | 9606.0              | CD300H    | 1.0               | 1.0              | MTQRAGAAMLPSALLLLCVPGCLTVSGPSTVMGAVGESLSVQCRYEEKYKTFNKYWCRQP |
