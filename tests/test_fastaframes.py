from dataclasses import asdict
from io import StringIO

import pandas as pd

from fastaframes import fasta_to_entries, FastaEntry, to_df, df_to_entries, to_fasta
from fastaframes.util import convert_to_best_datatype


def test_fasta_to_entries():
    fasta_content = \
        '>sp|A0A087X1C5|CP2D7_HUMAN Putative cytochrome P450 2D7 OS=Homo sapiens OX=9606 GN=CYP2D7 PE=5 SV=1' \
        '\nMGLEALVPLAMIVAIFLLLVDLMHR\nHQRWAARYPPGPLPLPGLGNLLH\nVDFQNTPYCFDQ\n'
    file_input = StringIO(fasta_content)
    result = fasta_to_entries(file_input)
    expected = [FastaEntry(
        db='sp',
        unique_identifier='A0A087X1C5',
        entry_name='CP2D7_HUMAN',
        protein_name='Putative cytochrome P450 2D7',
        organism_name='Homo sapiens',
        organism_identifier='9606',
        gene_name='CYP2D7',
        protein_existence='5',
        sequence_version='1',
        protein_sequence='MGLEALVPLAMIVAIFLLLVDLMHRHQRWAARYPPGPLPLPGLGNLLHVDFQNTPYCFDQ'
    )]
    assert result == expected


def test_fasta_to_entries_no_newlines():
    fasta_content = \
        '>sp|A0A087X1C5|CP2D7_HUMAN Putative cytochrome P450 2D7 OS=Homo sapiens OX=9606 GN=CYP2D7 PE=5 SV=1' \
        '\nMGLEALVPLAMIVAIFLLLVDLMHRHQRWAARYPPGPLPLPGLGNLLHVDFQNTPYCFDQ\n'
    file_input = StringIO(fasta_content)
    result = fasta_to_entries(file_input)
    expected = [FastaEntry(
        db='sp',
        unique_identifier='A0A087X1C5',
        entry_name='CP2D7_HUMAN',
        protein_name='Putative cytochrome P450 2D7',
        organism_name='Homo sapiens',
        organism_identifier='9606',
        gene_name='CYP2D7',
        protein_existence='5',
        sequence_version='1',
        protein_sequence='MGLEALVPLAMIVAIFLLLVDLMHRHQRWAARYPPGPLPLPGLGNLLHVDFQNTPYCFDQ'
    )]
    assert result == expected


def test_from_fasta():
    fasta_content = \
        '>sp|A0A087X1C5|CP2D7_HUMAN Putative cytochrome P450 2D7 OS=Homo sapiens OX=9606 GN=CYP2D7 PE=5 SV=1' \
        '\nMGLEALVPLAMIVAIFLLLVDLMHRHQRWAARYPPGPLPLPGLGNLLHVDFQNTPYCFDQ\n'
    file_input = StringIO(fasta_content)
    result = to_df(file_input)
    expected = pd.DataFrame([asdict(FastaEntry(
        db='sp',
        unique_identifier='A0A087X1C5',
        entry_name='CP2D7_HUMAN',
        protein_name='Putative cytochrome P450 2D7',
        organism_name='Homo sapiens',
        organism_identifier='9606',
        gene_name='CYP2D7',
        protein_existence='5',
        sequence_version='1',
        protein_sequence='MGLEALVPLAMIVAIFLLLVDLMHRHQRWAARYPPGPLPLPGLGNLLHVDFQNTPYCFDQ'
    ))])

    for col_name in expected:
        expected[col_name] = convert_to_best_datatype(expected[col_name])

    pd.testing.assert_frame_equal(result, expected)


def test_df_to_entries():
    data = {
        'db': ['sp'],
        'unique_identifier': ['A0A087X1C5'],
        'entry_name': ['CP2D7_HUMAN'],
        'protein_name': ['Putative cytochrome P450 2D7'],
        'organism_name': ['Homo sapiens'],
        'organism_identifier': ['9606'],
        'gene_name': ['CYP2D7'],
        'protein_existence': ['5'],
        'sequence_version': ['1'],
        'protein_sequence': ['MGLEALVPLAMIVAIFLLLVDLMHRHQRWAARYPPGPLPLPGLGNLLHVDFQNTPYCFDQ']
    }

    fasta_df = pd.DataFrame(data)
    result = df_to_entries(fasta_df)
    expected = [FastaEntry(
        db='sp',
        unique_identifier='A0A087X1C5',
        entry_name='CP2D7_HUMAN',
        protein_name='Putative cytochrome P450 2D7',
        organism_name='Homo sapiens',
        organism_identifier='9606',
        gene_name='CYP2D7',
        protein_existence='5',
        sequence_version='1',
        protein_sequence='MGLEALVPLAMIVAIFLLLVDLMHRHQRWAARYPPGPLPLPGLGNLLHVDFQNTPYCFDQ'
    )]
    assert result == expected


def test_to_fasta():
    data = {
        'db': ['sp'],
        'unique_identifier': ['A0A087X1C5'],
        'entry_name': ['CP2D7_HUMAN'],
        'protein_name': ['Putative cytochrome P450 2D7'],
        'organism_name': ['Homo sapiens'],
        'organism_identifier': ['9606'],
        'gene_name': ['CYP2D7'],
        'protein_existence': ['5'],
        'sequence_version': ['1'],
        'protein_sequence': ['MGLEALVPLAMIVAIFLLLVDLMHRHQRWAARYPPGPLPLPGLGNLLHVDFQNTPYCFDQ']
    }
    fasta_df = pd.DataFrame(data)
    result = to_fasta(fasta_df)
    expected = \
        '>sp|A0A087X1C5|CP2D7_HUMAN PN=Putative cytochrome P450 2D7 OS=Homo sapiens OX=9606 GN=CYP2D7 PE=5 SV=1' \
        '\nMGLEALVPLAMIVAIFLLLVDLMHRHQRWAARYPPGPLPLPGLGNLLHVDFQNTPYCFDQ\n'
    assert result.getvalue() == expected


def test_to_fasta_from_file():
    with open('tests/data/test.fasta', 'r') as file_input:
        result = fasta_to_entries(file_input)

    expected = [
        FastaEntry(
            db='sp',
            unique_identifier='A0A087X1C5',
            entry_name='CP2D7_HUMAN',
            protein_name='Putative cytochrome P450 2D7',
            organism_name='Homo sapiens',
            organism_identifier='9606',
            gene_name='CYP2D7',
            protein_existence='5',
            sequence_version='1',
            protein_sequence='MGLEALVPLAMIVAIFLLLVDLMHRHQRWAARYPPGPLPLPGLGNLLHVDFQNTPYCFDQ'),

        FastaEntry(
            db='sp',
            unique_identifier='A0A0B4J2F2',
            entry_name='SIK1B_HUMAN',
            protein_name='Putative serine/threonine-protein kinase SIK1B',
            organism_name='Homo sapiens',
            organism_identifier='9606',
            gene_name='SIK1B',
            protein_existence='5',
            sequence_version='1',
            protein_sequence='MVIMSEFSADPAGQGQGQQKPLRVGFYDIERTLGKGNFAVVKLARHRVTKTQVAIKIIDKLVQ'
        ),
        FastaEntry(
            db='sp',
            unique_identifier='A0A0C5B5G6',
            entry_name='MOTSC_HUMAN',
            protein_name='Mitochondrial-derived peptide MOTS-c',
            organism_name='Homo sapiens',
            organism_identifier='9606',
            gene_name='MT-RNR1',
            protein_existence='1',
            sequence_version='1',
            protein_sequence='MRWQEMGYIFYPRKLR'
        ),
        FastaEntry(
            db='sp',
            unique_identifier='A0A0K2S4Q6',
            entry_name='CD3CH_HUMAN',
            protein_name='Protein CD300H',
            organism_name='Homo sapiens',
            organism_identifier='9606',
            gene_name='CD300H',
            protein_existence='1',
            sequence_version='1',
            protein_sequence='MTQRAGAAMLPSALLLLCVPGCLTVSGPSTVMGAVGESLSVQCRYEEKYKTFNKYWCRQP'
        )
    ]

    assert result == expected

    io = to_fasta(result)
    result2 = fasta_to_entries(io)

    assert result2 == expected
