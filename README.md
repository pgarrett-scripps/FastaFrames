![example workflow](https://github.com/pgarrett-scripps/FastaFrames/actions/workflows/python-package.yml/badge.svg)
![example workflow](https://github.com/pgarrett-scripps/FastaFrames/actions/workflows/pylint.yml/badge.svg)

# FastaFrames
FastaFrames is a python package to convert between FASTA files and pandas DataFrames.

## Usage

To install fastaframes use pip:

```sh
pip install fastaframes
```

### Reading a FASTA file
```python
from fastaframes import to_df

fasta_df = to_df(data='example.fasta')
```

### Writing a FASTA file
```python
from fastaframes import to_fasta

to_fasta(data=fasta_df, output_file='output.fasta')
```

# Columns:
- **db**: Database from which the sequence was retrieved. db is 'sp' for UniProtKB/Swiss-Prot and 'tr' for UniProtKB/TrEMBL.
- **unique_identifier**: The primary accession number of the UniProtKB entry.
- **entry_name**: The entry name of the UniProtKB entry.
- **protein_name**: The recommended name of the UniProtKB entry as annotated in the RecName field. For UniProtKB/TrEMBL entries without a RecName field, the SubName field is used. In case of multiple SubNames, the first one is used. The 'precursor' attribute is excluded, 'Fragment' is included with the name if applicable.
- **organism_name**:  The scientific name of the organism of the UniProtKB entry.
- **organism_identifier**: The unique identifier of the source organism, assigned by the NCBI.
- **gene_name**: The first gene name of the UniProtKB entry. If there is no gene name, OrderedLocusName or ORFname, the GN field is not listed.
- **protein_existence**: The numerical value describing the evidence for the existence of the protein.
- **sequence_version**: The version number of the sequence.
- **protein_sequence**: The protein amino acid sequence.

## Example FASTA file:

```
>sp|A0A087X1C5|CP2D7_HUMAN Putative cytochrome P450 2D7 OS=Homo sapiens OX=9606 GN=CYP2D7 PE=5 SV=1
MGLEALVPLAMIVAIFLLLVDLMHRHQRWAARYPPGPLPLPGLGNLLHVDFQNTPYCFDQ
```

## Will produce the following:

|   | db | unique_identifier | entry_name   | protein_name                                         | organism_name | organism_identifier | gene_name | protein_existence | sequence_version | protein_sequence                                       |
|---|----|------------------|--------------|------------------------------------------------------|---------------|---------------------|-----------|-------------------|------------------|--------------------------------------------------------|
| 0 | sp | A0A087X1C5       | CP2D7_HUMAN  | Putative cytochrome P450 2D7                         | Homo sapiens  | 9606.0              | CYP2D7    | 5.0               | 1.0              | MGLEALVPLAMIVAIFLLLVDLMHRHQRWAARYPPGPLPLPGLGNLLHVDFQNTPYCFDQ |
