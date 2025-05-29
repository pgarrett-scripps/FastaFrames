import streamlit as st  # pip install streamlit
from src.fastaframes.fastaframes import to_df

fasta = st.file_uploader("Upload FASTA file", type=".fasta")


def get_protein_id(row):

    elems = [row['db'], row['unique_identifier'], row['entry_name']]
    elems = [e for e in elems if e is not None and e != 'None']
    return "|".join(elems)

if fasta:
    df = to_df(fasta)

    # 
    df['protein_id'] = df.apply(get_protein_id, axis=1)

    # display the dataframe
    st.dataframe(df)

    