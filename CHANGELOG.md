# Changelog

All notable changes to this project will be documented in this file.

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