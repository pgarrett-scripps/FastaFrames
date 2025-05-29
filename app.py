import streamlit as st  # pip install streamlit
from src.fastaframes.fastaframes import to_df
import requests

st.set_page_config(
    page_title="FASTA to DataFrame",
    page_icon=":dna:",
)

st.title("FASTA to DataFrame Converter")
example_fasta = "https://github.com/pgarrett-scripps/FastaFrames/blob/main/example.fasta"
st.caption("Upload a FASTA file and convert it to a Pandas DataFrame with structured metadata. FastaFrames parses " \
"UniProt-formatted FASTA files and extracts fields like database source, identifiers, protein names, organism details, " \
"gene names, and protein sequences into organized columns. Try this example: " \
f"[example.fasta]({example_fasta})")

fasta = st.file_uploader("Upload FASTA file", type=".fasta")

if fasta:
    df = to_df(fasta)

    # display the dataframe
    st.dataframe(df, use_container_width=True, hide_index=True)
    
    # Show download option
    st.download_button(
        label="Download as CSV",
        data=df.to_csv(index=False).encode('utf-8'),
        file_name=f'{fasta.name}.csv',
        on_click='ignore',
        use_container_width=True,
        help="Download the DataFrame as a CSV file.",
        mime='text/csv',
    )
else:
    st.info("Please upload a FASTA file")

