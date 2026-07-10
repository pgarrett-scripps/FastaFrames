# CLAUDE.md

Guidance for LLM coding agents working in this repository.

## What this package is

FastaFrames converts between UniProt FASTA files and pandas DataFrames (with the
`FastaEntry` dataclass in between). Keep it **simple and focused** — a "very
simple fasta file parser," not a general bioinformatics toolkit.

## Layout

```
src/fastaframes/
  __init__.py      # public API re-exports + __version__ (source of truth for version)
  fastaframes.py   # FastaEntry dataclass + all parse/serialize/convert functions
  util.py          # get_lines input normalization
  exceptions.py    # FastaFramesError hierarchy (FastaFormatError, InvalidInputError)
  py.typed         # PEP 561 marker — keep it; the package ships type hints
tests/             # pytest suite
examples/          # runnable recipes (decoy generation, streamlit app, tutorial)
```

## Tooling (matches CI — `.github/workflows/ci.yml`)

- **ruff** for lint and format: `ruff check src tests`, `ruff format --check src tests`.
- **ty** for type checking: `ty check src`.
- **pytest** on Python 3.10–3.13.
- `just check` runs lint + typecheck + test; `just fmt` formats. Dev install:
  `pip install -e ".[dev]"`.

## Conventions

- Target **Python >= 3.10**: use `str | None`, `list[...]`, `dict[...]`,
  `collections.abc` generics — NOT `typing.Optional`/`List`/`Dict`.
- ruff line length is 120; keep `__all__` isort-sorted (RUF022).
- Full Sphinx-style (`:param:`/`:return:`) docstrings on public functions.
- Bump `__version__` in `__init__.py` and add a `CHANGELOG.md` entry for any
  user-visible change (packaging reads the version from `__init__`).
- **Do not commit `.venv/`** or build artifacts (gitignored).

## Invariants — do not break these (each has a test)

- **Serialized output ends every record with a trailing `\n`** (`serialize`
  returns `header\nsequence\n`); `entries_to_fasta` concatenates records
  directly. `test_to_fasta` pins this.
- **Round-trip fidelity**, including `additional_fields` (non-standard
  `KEY=value` header fields), which are preserved through parse → serialize.
- **`header` / `protein_id` never emit the literal string `"None"`** — missing
  components render empty.
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
