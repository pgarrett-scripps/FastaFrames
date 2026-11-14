# Contributing to FastaFrames

Bug reports, fixes and new features are welcome. Use the
[issue tracker](https://github.com/pgarrett-scripps/FastaFrames/issues) for bugs and
questions, with a minimal example and the FastaFrames, fastatacular, pandas and Python
versions. Report security problems privately as described in [SECURITY.md](SECURITY.md).

FastaFrames is deliberately small: it reads FASTA files into pandas DataFrames and writes
them back. All FASTA parsing and writing is done by
[fastatacular](https://pypi.org/project/fastatacular/), and FastaFrames only maps its
`SequenceEntry` fields to DataFrame columns. A parser, a header format or a sequence
operation belongs in fastatacular, not here. If you are unsure which repository your change
belongs in, open an issue and ask before writing it.

## Claiming an issue

Comment on the issue to say you are taking it before you start, and check that nobody has
already done so. One issue, one contributor: when two pull requests arrive for the same issue,
the first one that was claimed is reviewed and the other is closed, however good it is. If a
claim has had no pull request or update for 14 days, comment to take it over.

Issues labelled `good first issue` need no knowledge of the codebase. Issues labelled
`help wanted` do: each one names the files to read first and the test that shows it is done.

## Development setup

You need Python 3.12 or newer and [just](https://just.systems/). FastaFrames requires 3.12
because fastatacular does.

```bash
git clone https://github.com/pgarrett-scripps/FastaFrames.git
cd FastaFrames
just install
```

`just install` runs `pip install -e ".[dev]"`, which installs the package in editable mode
plus the `dev` extra that the checks need (pytest, ruff and ty). Use a virtual environment.
FastaFrames is a standalone repository and is not part of the tacular-omics uv workspace, so
there is no `uv sync` step and no workspace root to run commands from: every command below is
run in this repository.

A correct install runs the whole suite with nothing skipped. Read the passed count at the end
of the run, not the collected total.

`just --list` shows the other recipes.

## Running checks

CI runs these on every pull request. Run them before you push:

```bash
just lint        # ruff check src tests
just fmt-check   # ruff format --check src tests
just typecheck   # ty check src
just test        # python -m pytest tests
```

`just check` runs all four in that order, and is the single command CI's jobs add up to.
`just fmt` rewrites the files that `just fmt-check` complains about, and `just fix` applies
ruff's automatic lint fixes.

Tests run on Python 3.12, 3.13 and 3.14 in CI. The suite is fast, under a second, so there is
no reason to push without having run it.

## Pull requests

- Branch from `main` and keep each pull request to one change.
- Add tests for new behaviour and fixed bugs. Prefer round-trip and multi-record assertions
  over single-value spot checks: past bugs in this package slipped through because the tests
  only covered the single-entry happy path.
- Add a line for user-visible changes at the top of `CHANGELOG.md`, under an
  `## [Unreleased]` heading (add the heading if it is not there yet).
- Do not change the version number. It lives in `src/fastaframes/__init__.py` and the
  maintainers set it when releasing.
- Do not change the `fastatacular` or `pandas` version specifiers in `pyproject.toml`. Open an
  issue instead: a floor or a cap change affects what every user gets and is a maintainer
  decision.
- CI must pass: lint, format check, type check, and tests on Python 3.12 to 3.14.

## Things not to break

Each of these has a test, and `CLAUDE.md` lists them with the test that pins each one:

- Every written line ends with `\n`, and an empty sequence writes its header line alone.
- UniProt files round-trip byte-identically through `to_df` then `to_fasta`, including
  `additional_fields`, NaN gaps and CSV round trips.
- No invented header text: no `PN=`, no `nan`, no `None`. An identifier that is not three
  `|`-separated fields stays whole in `unique_identifier` and is written back as it was read.
- Unknown header fields are kept in `additional_fields` rather than raising.
- All package errors derive from `FastaFramesError`. `FastaFormatError` stays a `ValueError`
  subclass and `InvalidInputError` stays a `TypeError`.
- Library code logs under `logging.getLogger(__name__)` and installs only a `NullHandler`.

## CLAUDE.md and AI agents

`CLAUDE.md` is the maintained developer guide: layout, tooling, conventions and the
invariants above. Read it before a larger change, and update it in the same pull request when
your change makes part of it wrong. AI-assisted contributions are welcome, and you are
responsible for reviewing and testing what you submit.
