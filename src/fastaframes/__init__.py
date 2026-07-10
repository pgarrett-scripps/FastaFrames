"""fastaframes - A simple fasta file parser.

Errors derive from :class:`FastaFramesError` (see :mod:`fastaframes.exceptions`).
The package logs under the ``"fastaframes"`` logger and installs a
:class:`logging.NullHandler`, so it stays silent until logging is configured::

    import logging
    logging.getLogger("fastaframes").setLevel(logging.DEBUG)
"""

import logging

from .exceptions import (
    FastaFormatError,
    FastaFramesError,
    InvalidInputError,
)
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

# Library best practice: a NullHandler so importing fastaframes never emits
# logging output unless the application explicitly configures handlers.
logging.getLogger(__name__).addHandler(logging.NullHandler())

__all__ = [
    "FastaEntry",
    "FastaFields",
    "FastaFormatError",
    "FastaFramesError",
    "InvalidInputError",
    "df_to_entries",
    "entries_to_df",
    "entries_to_fasta",
    "fasta_to_entries",
    "to_df",
    "to_fasta",
]

__version__ = "1.4.0"
