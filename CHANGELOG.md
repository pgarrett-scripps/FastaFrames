# Changelog

All notable changes to this project will be documented in this file.

## [1.4.0]

### Added
- `FastaEntry.additional_fields`: non-standard `KEY=value` header fields (beyond
  PN/OS/OX/GN/PE/SV) are now preserved and round-tripped instead of raising.
- Optional sequence line wrapping via the `line_wrapping` and
  `max_sequence_length` arguments to `serialize()` / `entries_to_fasta()` /
  `to_fasta()`. `max_sequence_length=None`/`0` means no wrapping.
- Exception hierarchy in `fastaframes.exceptions`, exported from the package
  root: `FastaFramesError` (base), `FastaFormatError` (a `ValueError`, with
  `.header`/`.reason`), and `InvalidInputError` (a `TypeError`).
- Structured logging under the `fastaframes` logger (with a `NullHandler`).
- Expanded docstrings across the public API.

### Fixed
- `FastaEntry.header` renders missing `db`/`unique_identifier`/`entry_name` as
  empty strings instead of the literal text `None`.
- `df_to_entries` tolerates missing columns (filling defaults with a warning)
  instead of raising `KeyError`.

## [1.3.0]
- Production-readiness overhaul: modern packaging (`pyproject.toml`, removed
  `setup.py`), `ruff` + `ty` tooling, consolidated CI, Python 3.10+ typing, and
  a rewritten README.

## [1.2.2]
- fixed None being cast to 'None' in dataframe
- added protein_id col to df
- removed convert_to_best_datatype

## [1.2.1]
- bugs
- streamlit community cloud app
- black formatting
- updated reqs with pipreqs

## [1.2.1]
- fix for malformed fastas (warning instead of error)

## [1.0.0]

## Added
- github open source project requirements
- examples

## Changed
- improved readability of get_lines
- added Enum for FASTA info

## [0.0.3]

## Changes
- df_to_entries now filters by expected columns prior to creating entries
- _get_lines now works with streamlit uploaded file, and any io-type 

## [0.0.2]

## Added
- to_df and to_fasta: versatile functions 
- example.py
- example.fasta

## Changes
- Better documentation
- Simpler README
- Moved FastaEntry serialize function within class

## [0.0.1]

### Added
- Functions to handle basic parsing to/from fasta file, FastaEntry dataclasses, and pandas dataframes
- Test Suite