"""
Exception hierarchy for fastaframes.

Every package-specific failure derives from :class:`FastaFramesError`, so a
caller (including an automated agent) can catch them all with a single
``except FastaFramesError``. Each concrete exception also derives from the
closest built-in (``ValueError`` / ``TypeError``) so code that catches those
broad types keeps working, and carries structured attributes so callers can
branch on *what* went wrong without matching message text.
"""


class FastaFramesError(Exception):
    """Base class for every error raised by fastaframes."""


class FastaFormatError(FastaFramesError, ValueError):
    """A FASTA header or record could not be parsed.

    Also a :class:`ValueError` for backward compatibility.

    :ivar header: The offending header/record text.
    :ivar reason: Human-readable explanation of the failure.
    """

    def __init__(self, header: object, reason: str = "could not be parsed as a FASTA header") -> None:
        self.header = header
        self.reason = reason
        # The "Invalid fasta header format" prefix is part of the public
        # contract (callers/tests match on it) - keep it stable.
        super().__init__(f"Invalid fasta header format: {header} ({reason})")


class InvalidInputError(FastaFramesError, TypeError):
    """An input value was not a supported source of FASTA lines.

    Raised by :func:`fastaframes.util.get_lines` for inputs that are neither a
    path/string, a readable file-like object, nor an iterable of lines. Also a
    :class:`TypeError`.

    :ivar value: The unsupported input value.
    :ivar received_type: The ``type`` of the unsupported input value.
    """

    def __init__(self, value: object) -> None:
        self.value = value
        self.received_type = type(value)
        super().__init__(
            f"Unsupported input type {self.received_type.__name__!r}: expected a file path, a FASTA string, "
            "a readable file-like object, or an iterable of str/bytes lines."
        )
