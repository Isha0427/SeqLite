import streamlit as st

st.set_page_config(page_title="Bioinformatics Toolkit", layout="wide")
st.title("🧬 SeqLite ")

# Sidebar Navigation
tools = [
    "Home",
    "Sequence Translator",
    "GC Content Calculator",
    "Reverse Complement",
    "Protein Molecular Weight",
    "About"
]

selection = st.sidebar.radio("Select a Tool", tools)

# -------------------- HOME --------------------
if selection == "Home":
    st.header("Welcome to the Bioinformatics Toolkit")
    st.markdown("""
    This toolkit provides a collection of commonly used bioinformatics utilities for students and researchers:

    - DNA/RNA Sequence tools
    - Protein analysis
    - Interactive results

    Select a tool from the sidebar to begin.
    """)

# -------------------- SEQUENCE TRANSLATOR --------------------
elif selection == "Sequence Translator":
    st.header("🧬 Sequence Translator")
    dna_seq = st.text_area("Enter DNA Sequence:")

    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

    if dna_seq:
        dna_seq = dna_seq.upper().replace(" ", "").replace("\n", "")
        protein_seq = ""
        for i in range(0, len(dna_seq) - 2, 3):
            codon = dna_seq[i:i+3]
            protein_seq += codon_table.get(codon, 'X')
        st.text_area("Protein Sequence:", value=protein_seq, height=150)

# -------------------- GC CONTENT CALCULATOR --------------------
elif selection == "GC Content Calculator":
    st.header("🧪 GC Content Calculator")
    sequence = st.text_area("Enter DNA/RNA Sequence:")
    if sequence:
        sequence = sequence.upper()
        g = sequence.count("G")
        c = sequence.count("C")
        gc_content = (g + c) / len(sequence) * 100 if len(sequence) > 0 else 0
        st.write(f"GC Content: **{gc_content:.2f}%**")

# -------------------- REVERSE COMPLEMENT --------------------
elif selection == "Reverse Complement":
    st.header("🔄 Reverse Complement Tool")
    dna = st.text_area("Enter DNA Sequence:")
    comp_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    if dna:
        dna = dna.upper()
        try:
            rev_comp = ''.join([comp_dict[base] for base in reversed(dna)])
            st.text_area("Reverse Complement:", value=rev_comp, height=150)
        except KeyError:
            st.error("Invalid nucleotide detected. Use only A, T, G, C.")

# -------------------- PROTEIN MOLECULAR WEIGHT --------------------
elif selection == "Protein Molecular Weight":
    st.header("⚖️ Protein Molecular Weight Calculator")
    aa_weights = {
        'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10,
        'C': 121.16, 'Q': 146.15, 'E': 147.13, 'G': 75.07,
        'H': 155.16, 'I': 131.18, 'L': 131.18, 'K': 146.19,
        'M': 149.21, 'F': 165.19, 'P': 115.13, 'S': 105.09,
        'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
    }
    protein = st.text_area("Enter Protein Sequence:")
    if protein:
        protein = protein.upper().strip()
        try:
            mw = sum(aa_weights[aa] for aa in protein)
            st.write(f"Molecular Weight: **{mw:.2f} Da**")
        except KeyError:
            st.error("Invalid amino acid detected. Use 1-letter codes only (A, R, N, etc.)")

# -------------------- ABOUT --------------------
elif selection == "About":
   st.header("📘 About Page")
   st.markdown("""
    ### 🔬 Overview
    This webpage provides information of DNA sequence .

    ### 🎯 Objectives
    - Help users understand DNA/RNA sequence
    - Provide quick access to GC content , Protein Molecular Weight ,etc
    - Easy input and fast results

    ### 🚀 Features
    - Clean user interface
    - Instant response
    """)

   st.divider()
   st.header("👥 Team Page")
   st.markdown("""
    ### 👩‍💻 Contributor
    - **Name:** Isha Jitendra Hajare
    - **Role:** Bioinformatics Student
    - **Contact:** ishahajare27@gmail.com
    """)
