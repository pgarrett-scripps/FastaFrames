# Changelog

All notable changes to this project will be documented in this file.

## [2.0.0]

FastaFrames is now a pandas layer over
[fastatacular](https://pypi.org/project/fastatacular/), which does all FASTA
parsing and writing. The public API and DataFrame columns are unchanged.

### Breaking
- Requires Python >= 3.12 (fastatacular's floor), `fastatacular>=1.3,<2` and
  `pandas>=2.1`.
- Written headers changed. The protein name is no longer written as an invented
  `PN=` field, and a header whose identifier is not `db|accession|entry_name`
  is written as it was read: `>myid some desc` stays `>myid some desc` instead
  of becoming `>|myid| PN=some desc`. `FastaEntry.header` for an entry without
  `db` and `entry_name` is `>ACCESSION`, not `>|ACCESSION|`. Files written by
  1.x still read correctly: a leading `PN=` is read into `protein_name`.
- A `str` that is not FASTA text (no newline, does not start with `>`) is a
  path; a missing path raises `FileNotFoundError` instead of returning an
  empty DataFrame.
- Sequence lines before the first header raise `FastaFormatError` (skipped with
  a warning under `skip_error=True`) instead of being dropped silently.
- Whitespace around the protein name and field values is trimmed.
- A header that is not three `|`-separated fields no longer emits a
  `UserWarning`; it is kept whole in `unique_identifier` as before.
- A record with an empty sequence is written as its header line alone.
- Writing an entry with no identifier, or a non-integer `organism_identifier`,
  `protein_existence` or `sequence_version`, raises `fastatacular.FastaWriteError`
  (a `ValueError`); every entry is checked before anything is written.

### Changed
- Header parsing and FASTA writing go through fastatacular; FastaFrames' own
  header parser is removed.
- Paths may be `str` or `pathlib.Path`, plain or gzip/bzip2/xz compressed;
  `output_file` ending in `.gz`, `.bz2` or `.xz` is written compressed.

### Fixed
- UniProt files round-trip byte-identically through `to_df` -> `to_fasta`.
- Missing values (NaN, None, `pd.NA`) are left out of written headers instead of
  being written as `nan` (`GN=nan`, `>nan|myid|nan`); integral floats from a
  column with gaps or a CSV round trip (`9606.0`) are written as `9606`.
- `additional_fields` survive a DataFrame round trip (also through CSV, where
  the column comes back as text).
- A `>` inside a header is kept (it was deleted).
- A tab after the identifier separates it from the description.
- A `.gz` path no longer crashes; `pathlib.Path` inputs are accepted.
- A NaN sequence with `max_sequence_length` set no longer raises `TypeError`.

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