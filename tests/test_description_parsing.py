from dataclasses import asdict

import pytest
from fastaframes.fastaframes import _fasta_str_to_entry, fasta_to_entries, FastaEntry


@pytest.mark.parametrize("fasta_entry, expected", [
    (">sp|Q8I6R7|ACN2_ACAGO Acanthoscurrin-2 (Fragment) OS=Acanthoscurria gomesiana OX=115339 GN=acantho2 PE=1 SV=1", {
        'db': 'sp',
        'unique_identifier': 'Q8I6R7',
        'entry_name': 'ACN2_ACAGO',
        'protein_name': 'Acanthoscurrin-2 (Fragment)',
        'organism_name': 'Acanthoscurria gomesiana',
        'organism_identifier': '115339',
        'gene_name': 'acantho2',
        'protein_existence': '1',
        'sequence_version': '1',
        'protein_sequence': '',
    }),
    (
    ">sp|P27748|ACOX_CUPNH Acetoin catabolism protein X OS=Cupriavidus necator (strain ATCC 17699 / H16 / DSM 428 / Stanier 337) OX=381666 GN=acoX PE=4 SV=2",
    {
        'db': 'sp',
        'unique_identifier': 'P27748',
        'entry_name': 'ACOX_CUPNH',
        'protein_name': 'Acetoin catabolism protein X',
        'organism_name': 'Cupriavidus necator (strain ATCC 17699 / H16 / DSM 428 / Stanier 337)',
        'organism_identifier': '381666',
        'gene_name': 'acoX',
        'protein_existence': '4',
        'sequence_version': '2',
        'protein_sequence': '',
    }),
    (">tr|Q3SA23|Q3SA23_9HIV1 Protein Nef (Fragment) OS=Human immunodeficiency virus 1  OX=11676 GN=nef PE=3 SV=1", {
        'db': 'tr',
        'unique_identifier': 'Q3SA23',
        'entry_name': 'Q3SA23_9HIV1',
        'protein_name': 'Protein Nef (Fragment)',
        'organism_name': 'Human immunodeficiency virus 1 ',
        'organism_identifier': '11676',
        'gene_name': 'nef',
        'protein_existence': '3',
        'sequence_version': '1',
        'protein_sequence': '',
    }),
    (">sp|A0A087X1C5|CP2D7_HUMAN Putative cytochrome P450 2D7 OS=Homo sapiens OX=9606 GN=CYP2D7 PE=5 SV=1", {
        'db': 'sp',
        'unique_identifier': 'A0A087X1C5',
        'entry_name': 'CP2D7_HUMAN',
        'protein_name': 'Putative cytochrome P450 2D7',
        'organism_name': 'Homo sapiens',
        'organism_identifier': '9606',
        'gene_name': 'CYP2D7',
        'protein_existence': '5',
        'sequence_version': '1',
        'protein_sequence': '',
    }),
    (">sp|A0A087X1C5|CP2D7_HUMAN Putative cytochrome P450 2D7 OS=Homo sapiens OX=9606 PE=5 SV=1", {
        'db': 'sp',
        'unique_identifier': 'A0A087X1C5',
        'entry_name': 'CP2D7_HUMAN',
        'protein_name': 'Putative cytochrome P450 2D7',
        'organism_name': 'Homo sapiens',
        'organism_identifier': '9606',
        'gene_name': None,
        'protein_existence': '5',
        'sequence_version': '1',
        'protein_sequence': '',
    }),
    (">tr|G3MXS6|G3MXS6_BOVIN Uncharacterized protein (Fragment) OS=Bos taurus PE=4 SV=1", {
        'db': 'tr',
        'unique_identifier': 'G3MXS6',
        'entry_name': 'G3MXS6_BOVIN',
        'protein_name': 'Uncharacterized protein (Fragment)',
        'organism_name': 'Bos taurus',
        'organism_identifier': None,
        'gene_name': None,
        'protein_existence': '4',
        'sequence_version': '1',
        'protein_sequence': '',
    }),
    (">tr|G3MXS6|G3MXS6_BOVIN Uncharacterized protein (Fragment) PE=4 SV=1", {
        'db': 'tr',
        'unique_identifier': 'G3MXS6',
        'entry_name': 'G3MXS6_BOVIN',
        'protein_name': 'Uncharacterized protein (Fragment)',
        'organism_name': None,
        'organism_identifier': None,
        'gene_name': None,
        'protein_existence': '4',
        'sequence_version': '1',
        'protein_sequence': '',
    }),
    (">tr|G3MXS6|G3MXS6_BOVIN Uncharacterized protein (Fragment)  PE=4 SV=1", { # note the extra space
        'db': 'tr',
        'unique_identifier': 'G3MXS6',
        'entry_name': 'G3MXS6_BOVIN',
        'protein_name': 'Uncharacterized protein (Fragment) ', # note the extra space
        'organism_name': None,
        'organism_identifier': None,
        'gene_name': None,
        'protein_existence': '4',
        'sequence_version': '1',
        'protein_sequence': '',
    }),
    (">tr|G3MXS6|G3MXS6_BOVIN", {
        'db': 'tr',
        'unique_identifier': 'G3MXS6',
        'entry_name': 'G3MXS6_BOVIN',
        'protein_name': None,
        'organism_name': None,
        'organism_identifier': None,
        'gene_name': None,
        'protein_existence': None,
        'sequence_version': None,
        'protein_sequence': '',
    })
,
    (">tr|G3MXS6|G3MXS6_BOVIN PE=4", {
        'db': 'tr',
        'unique_identifier': 'G3MXS6',
        'entry_name': 'G3MXS6_BOVIN',
        'protein_name': None,
        'organism_name': None,
        'organism_identifier': None,
        'gene_name': None,
        'protein_existence': '4',
        'sequence_version': None,
        'protein_sequence': '',
    })

    # Add more test cases if needed
])

def test_extract_fasta_info(fasta_entry, expected):
    result = asdict(_fasta_str_to_entry(fasta_entry))
    assert result == expected

