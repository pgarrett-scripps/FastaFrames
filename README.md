![example workflow](https://github.com/pgarrett-scripps/FastaFrames/actions/workflows/python-package.yml/badge.svg)
![example workflow](https://github.com/pgarrett-scripps/FastaFrames/actions/workflows/pylint.yml/badge.svg)

# FastaFrames
This Python module provides a set of functions to work with FASTA files. 
It allows you to read FASTA files, convert them to pandas dataframes, manipulate data, 
and write data back to FASTA files. It also supports converting FASTA files to a list of FastaEntry dataclass objects.

## Features
- Read FASTA files into pandas DataFrames
- Convert FASTA files to a list of FastaEntry objects
- Write FASTA files from pandas DataFrames
- Write FASTA files from a list of FastaEntry objects

## Usage
### Reading FASTA files
To read a FASTA file and convert it to a pandas DataFrame:

```python
from fastaframes import fasta_to_df

with open('example.fasta', 'r') as file_input:
    fasta_df = fasta_to_df(file_input)

print(fasta_df.head())
```

To read a FASTA file and convert it to a list of FastaEntry objects:

```python
from fastaframes import fasta_to_entries

with open('example.fasta', 'r') as file_input:
    entries = fasta_to_entries(file_input)

print(entries[:5])
```

### Writing FASTA files
To write a pandas DataFrame to a FASTA file:

```python
import pandas as pd
from fastaframes import to_fasta

fasta_df = pd.DataFrame() # empty
fasta_content = to_fasta(fasta_df)

with open('output.fasta', 'w') as output_file:
    output_file.write(fasta_content.getvalue())
```

To write a list of FastaEntry objects to a FASTA file:

```python
from fastaframes import to_fasta

entries = [] # empty
fasta_content = to_fasta(entries)

with open('output.fasta', 'w') as output_file:
    output_file.write(fasta_content.getvalue())
```

## Customization
The FastaEntry dataclass can be customized to store additional information or modify existing attributes
as needed. This can be done by editing the FastaEntry dataclass definition and updating the extract_fasta_info
function accordingly.