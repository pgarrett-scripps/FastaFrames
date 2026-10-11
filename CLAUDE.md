# CLAUDE.md

Guidance for LLM coding agents working in this repository.

## What this package is

FastaFrames converts between UniProt FASTA files and pandas DataFrames (with the
`FastaEntry` dataclass in between). Since 2.0 all FASTA parsing and writing is
done by **fastatacular** (`FastaReader`, `read_fasta`, `from_records`,
`write_fasta`); FastaFrames only maps `SequenceEntry` fields to its columns
(`_from_sequence_entry`, `_to_record`). Do not add a header parser here. Keep it **simple and focused** — a "very
simple fasta file parser," not a general bioinformatics toolkit.

## Layout

```
src/fastaframes/
  __init__.py      # public API re-exports + __version__ (source of truth for version)
  fastaframes.py   # FastaEntry dataclass + all parse/serialize/convert functions
  util.py          # as_source / get_lines input normalization (paths, text, handles)
  exceptions.py    # FastaFramesError hierarchy (FastaFormatError, InvalidInputError)
  py.typed         # PEP 561 marker — keep it; the package ships type hints
tests/             # pytest suite
examples/          # runnable recipes (decoy generation, streamlit app, tutorial)
```

## Tooling (matches CI — `.github/workflows/ci.yml`)

- **ruff** for lint and format: `ruff check src tests`, `ruff format --check src tests`.
- **ty** for type checking: `ty check src`.
- **pytest** on Python 3.12–3.14 (fastatacular requires 3.12).
- `just check` runs lint + typecheck + test; `just fmt` formats. Dev install:
  `pip install -e ".[dev]"`.

## Conventions

- Target **Python >= 3.12**: use `str | None`, `list[...]`, `dict[...]`,
  `collections.abc` generics — NOT `typing.Optional`/`List`/`Dict`.
- ruff line length is 120; keep `__all__` isort-sorted (RUF022).
- Full Sphinx-style (`:param:`/`:return:`) docstrings on public functions.
- Bump `__version__` in `__init__.py` and add a `CHANGELOG.md` entry for any
  user-visible change (packaging reads the version from `__init__`).
- **Do not commit `.venv/`** or build artifacts (gitignored).

## Invariants — do not break these (each has a test)

- **Every written line ends with `\n`**; an empty sequence writes the header
  line alone. `test_to_fasta` pins this.
- **UniProt files round-trip byte-identically** through `to_df` -> `to_fasta`,
  including `additional_fields`, NaN gaps and CSV round trips
  (`tests/test_regressions_2_0.py`).
- **No invented header text**: no `PN=`, no `nan`, no `None`. An identifier that
  is not three `|` fields stays whole in `unique_identifier` and is written back
  as is.
- **Line wrapping is parameter-only**: `serialize`/`entries_to_fasta`/`to_fasta`
  take `line_wrapping` + `max_sequence_length`; `None`/`0` means no wrap. There
  is intentionally **no `max_sequence_length` field** on `FastaEntry` (it would
  break the `asdict(...)` equality tests in `test_description_parsing`).
- **Unknown header fields are stored in `additional_fields`, not raised.**
- **`df_to_entries` tolerates missing columns** (warns, fills defaults).

## Errors and logging

- All package errors derive from `FastaFramesError` (`exceptions.py`).
  `FastaFormatError` must remain a `ValueError` subclass and keep the
  `"Invalid fasta header format"` message prefix. `InvalidInputError` is a
  `TypeError`. Raise `FastaFormatError` for malformed headers so
  `fasta_to_entries`' `skip_error` path catches the right type.
- Log under `logging.getLogger(__name__)` per module: normal flow → `DEBUG`,
  tolerated bad input → `WARNING`. The `NullHandler` lives in `__init__.py`;
  never add a non-Null handler in library code.

## When you change behavior

Add a test that exercises the actual failure mode. Prefer round-trip and
multi-record assertions over single-value spot checks — past bugs slipped
through because tests only covered the single-entry happy path.
