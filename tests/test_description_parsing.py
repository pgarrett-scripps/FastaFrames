from dataclasses import asdict

import pandas as pd
import pytest

from fastaframes.fastaframes import _fasta_str_to_entry, to_df


@pytest.mark.parametrize(
    "fasta_entry, expected",
    [
        (
            ">sp|Q8I6R7|ACN2_ACAGO Acanthoscurrin-2 (Fragment) OS=Acanthoscurria gomesiana OX=115339 GN=acantho2 PE=1 SV=1",
            {
                "db": "sp",
                "unique_identifier": "Q8I6R7",
                "entry_name": "ACN2_ACAGO",
                "protein_name": "Acanthoscurrin-2 (Fragment)",
                "organism_name": "Acanthoscurria gomesiana",
                "organism_identifier": "115339",
                "gene_name": "acantho2",
                "protein_existence": "1",
                "sequence_version": "1",
                "protein_sequence": "",
                "additional_fields": {},
            },
        ),
        (
            ">sp|P27748|ACOX_CUPNH Acetoin catabolism protein X OS=Cupriavidus necator (strain ATCC 17699 / H16 / DSM 428 / Stanier 337) OX=381666 GN=acoX PE=4 SV=2",
            {
                "db": "sp",
                "unique_identifier": "P27748",
                "entry_name": "ACOX_CUPNH",
                "protein_name": "Acetoin catabolism protein X",
                "organism_name": "Cupriavidus necator (strain ATCC 17699 / H16 / DSM 428 / Stanier 337)",
                "organism_identifier": "381666",
                "gene_name": "acoX",
                "protein_existence": "4",
                "sequence_version": "2",
                "protein_sequence": "",
                "additional_fields": {},
            },
        ),
        (
            ">tr|Q3SA23|Q3SA23_9HIV1 Protein Nef (Fragment) OS=Human immunodeficiency virus 1  OX=11676 GN=nef PE=3 SV=1",
            {
                "db": "tr",
                "unique_identifier": "Q3SA23",
                "entry_name": "Q3SA23_9HIV1",
                "protein_name": "Protein Nef (Fragment)",
                "organism_name": "Human immunodeficiency virus 1",
                "organism_identifier": "11676",
                "gene_name": "nef",
                "protein_existence": "3",
                "sequence_version": "1",
                "protein_sequence": "",
                "additional_fields": {},
            },
        ),
        (
            ">sp|A0A087X1C5|CP2D7_HUMAN Putative cytochrome P450 2D7 OS=Homo sapiens OX=9606 GN=CYP2D7 PE=5 SV=1",
            {
                "db": "sp",
                "unique_identifier": "A0A087X1C5",
                "entry_name": "CP2D7_HUMAN",
                "protein_name": "Putative cytochrome P450 2D7",
                "organism_name": "Homo sapiens",
                "organism_identifier": "9606",
                "gene_name": "CYP2D7",
                "protein_existence": "5",
                "sequence_version": "1",
                "protein_sequence": "",
                "additional_fields": {},
            },
        ),
        (
            ">sp|A0A087X1C5|CP2D7_HUMAN Putative cytochrome P450 2D7 OS=Homo sapiens OX=9606 PE=5 SV=1",
            {
                "db": "sp",
                "unique_identifier": "A0A087X1C5",
                "entry_name": "CP2D7_HUMAN",
                "protein_name": "Putative cytochrome P450 2D7",
                "organism_name": "Homo sapiens",
                "organism_identifier": "9606",
                "gene_name": None,
                "protein_existence": "5",
                "sequence_version": "1",
                "protein_sequence": "",
                "additional_fields": {},
            },
        ),
        (
            ">tr|G3MXS6|G3MXS6_BOVIN Uncharacterized protein (Fragment) OS=Bos taurus PE=4 SV=1",
            {
                "db": "tr",
                "unique_identifier": "G3MXS6",
                "entry_name": "G3MXS6_BOVIN",
                "protein_name": "Uncharacterized protein (Fragment)",
                "organism_name": "Bos taurus",
                "organism_identifier": None,
                "gene_name": None,
                "protein_existence": "4",
                "sequence_version": "1",
                "protein_sequence": "",
                "additional_fields": {},
            },
        ),
        (
            ">tr|G3MXS6|G3MXS6_BOVIN Uncharacterized protein (Fragment) PE=4 SV=1",
            {
                "db": "tr",
                "unique_identifier": "G3MXS6",
                "entry_name": "G3MXS6_BOVIN",
                "protein_name": "Uncharacterized protein (Fragment)",  # surrounding whitespace is trimmed
                "organism_name": None,
                "organism_identifier": None,
                "gene_name": None,
                "protein_existence": "4",
                "sequence_version": "1",
                "protein_sequence": "",
                "additional_fields": {},
            },
        ),
        (
            ">tr|G3MXS6|G3MXS6_BOVIN Uncharacterized protein (Fragment)  PE=4 SV=1",  # note the extra space
            {
                "db": "tr",
                "unique_identifier": "G3MXS6",
                "entry_name": "G3MXS6_BOVIN",
                "protein_name": "Uncharacterized protein (Fragment)",
                "organism_name": None,
                "organism_identifier": None,
                "gene_name": None,
                "protein_existence": "4",
                "sequence_version": "1",
                "protein_sequence": "",
                "additional_fields": {},
            },
        ),
        (
            ">tr|G3MXS6|G3MXS6_BOVIN",
            {
                "db": "tr",
                "unique_identifier": "G3MXS6",
                "entry_name": "G3MXS6_BOVIN",
                "protein_name": None,
                "organism_name": None,
                "organism_identifier": None,
                "gene_name": None,
                "protein_existence": None,
                "sequence_version": None,
                "protein_sequence": "",
                "additional_fields": {},
            },
        ),
        (
            ">tr|G3MXS6|G3MXS6_BOVIN PE=4",
            {
                "db": "tr",
                "unique_identifier": "G3MXS6",
                "entry_name": "G3MXS6_BOVIN",
                "protein_name": None,
                "organism_name": None,
                "organism_identifier": None,
                "gene_name": None,
                "protein_existence": "4",
                "sequence_version": None,
                "protein_sequence": "",
                "additional_fields": {},
            },
        ),
        (
            ">sp|Q8I6R7|ACN2_ACAGO Some protein OS=Homo sapiens XT=custom_value",
            {
                "db": "sp",
                "unique_identifier": "Q8I6R7",
                "entry_name": "ACN2_ACAGO",
                "protein_name": "Some protein",
                "organism_name": "Homo sapiens",
                "organism_identifier": None,
                "gene_name": None,
                "protein_existence": None,
                "sequence_version": None,
                "protein_sequence": "",
                "additional_fields": {"XT": "custom_value"},
            },
        ),
        (
            ">NP_000549.1 hemoglobin subunit alpha [Homo sapiens]",
            {
                "db": None,
                "unique_identifier": "NP_000549.1",
                "entry_name": None,
                "protein_name": "hemoglobin subunit alpha [Homo sapiens]",
                "organism_name": "Homo sapiens",
                "organism_identifier": None,
                "gene_name": None,
                "protein_existence": None,
                "sequence_version": None,
                "protein_sequence": "",
                "additional_fields": {},
            },
        ),
        (
            ">UniRef90_P69905 Hemoglobin subunit alpha n=12 Tax=Primates TaxID=9443 RepID=HBA_HUMAN",
            {
                "db": None,
                "unique_identifier": "UniRef90_P69905",
                "entry_name": None,
                "protein_name": "Hemoglobin subunit alpha",
                "organism_name": "Primates",
                "organism_identifier": "9443",
                "gene_name": None,
                "protein_existence": None,
                "sequence_version": None,
                "protein_sequence": "",
                "additional_fields": {
                    "n": "12",
                    "RepID": "HBA_HUMAN",
                },
            },
        ),
        (
            ">ENSZVIP00000000002.1 pep UG_Zviv_1:LG1:5289:19862:-1 gene:ENSZVIG00000000002.1 transcript:ENSZVIT00000000002.1 gene_biotype:protein_coding transcript_biotype:protein_coding",
            {
                "db": None,
                "unique_identifier": "ENSZVIP00000000002.1",
                "entry_name": None,
                "protein_name": "pep UG_Zviv_1:LG1:5289:19862:-1 gene:ENSZVIG00000000002.1 transcript:ENSZVIT00000000002.1 gene_biotype:protein_coding transcript_biotype:protein_coding",
                "organism_name": None,
                "organism_identifier": None,
                "gene_name": None,
                "protein_existence": None,
                "sequence_version": None,
                "protein_sequence": "",
                "additional_fields": {},
            },
        ),
    ],
)
def test_extract_fasta_info(fasta_entry: str, expected: dict[str, object]):
    result = asdict(_fasta_str_to_entry(fasta_entry))
    assert result.pop("raw_header") == fasta_entry[1:]
    assert result == expected



def test_to_df_non_uniprot_headers():
    fasta_content = """>NP_000549.1 hemoglobin subunit alpha [Homo sapiens]
MVLSPADKTN
>UniRef90_P69905 Hemoglobin subunit alpha n=12 Tax=Primates TaxID=9443 RepID=HBA_HUMAN
MVLSPADKTN
>ENSZVIP00000000002.1 pep UG_Zviv_1:LG1:5289:19862:-1 gene:ENSZVIG00000000002.1 transcript:ENSZVIT00000000002.1 gene_biotype:protein_coding transcript_biotype:protein_coding
MVLSPADKTN
"""

    df = to_df(fasta_content)

    assert {"organism_name", "organism_identifier"}.issubset(df.columns)
    assert len(df) == 3

    assert df.loc[0, "organism_name"] == "Homo sapiens"
    assert pd.isna(df.loc[0, "organism_identifier"])

    assert df.loc[1, "organism_name"] == "Primates"
    assert df.loc[1, "organism_identifier"] == "9443"

    assert pd.isna(df.loc[2, "organism_name"])
    assert pd.isna(df.loc[2, "organism_identifier"])
