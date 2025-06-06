import streamlit as st  # pip install streamlit
from fastaframes import to_df

fasta = st.file_uploader("Upload FASTA file", type=".fasta")

if fasta:
    df = to_df(fasta)

    # display the dataframe
    st.dataframe(df)

    st.balloons()  # 🎈
