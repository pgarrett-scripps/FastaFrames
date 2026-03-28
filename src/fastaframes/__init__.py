"""fastaframes - A simple fasta file parser."""

from .fastaframes import (
    FastaEntry,
    FastaFields,
    df_to_entries,
    entries_to_df,
    entries_to_fasta,
    fasta_to_entries,
    to_df,
    to_fasta,
)

__all__ = [
    "FastaEntry",
    "FastaFields",
    "df_to_entries",
    "entries_to_df",
    "entries_to_fasta",
    "fasta_to_entries",
    "to_df",
    "to_fasta",
]

__version__ = "1.3.0"
