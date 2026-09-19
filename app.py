import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt

from config import APP_TITLE, APP_ICON
from services.ncbi import get_human_mtco1, get_protein, get_homologs
from analysis.dna import nucleotide_stats, sliding_gc
from analysis.protein import protein_stats, amino_acid_counts, hydrophobicity
from analysis.alignment import pairwise_identity_matrix, global_alignment
from analysis.phylogeny import build_nj_tree
from analysis.variants import analyze_variant

st.set_page_config(page_title=APP_TITLE, page_icon=APP_ICON, layout="wide")

@st.cache_data(ttl=86400)
def load_data():
    human = get_human_mtco1()
    protein = get_protein()
    homologs = get_homologs()
    return human, protein, homologs

try:
    human, refseq_protein, homologs = load_data()
except Exception as e:
    st.error("NCBI data could not be loaded.")
    st.exception(e)
    st.stop()

dna = human["sequence"]
protein = refseq_protein["sequence"]
dna_stats = nucleotide_stats(dna)
prot_stats = protein_stats(protein)

st.sidebar.title("🧬 MitoGene Analyzer")
st.sidebar.caption("Real NCBI/RefSeq data • educational use")
page = st.sidebar.radio(
    "Navigation",
    [
        "🏠 Dashboard",
        "🧬 Gene Information",
        "🧪 DNA Analysis",
        "🧫 Protein Analysis",
        "🔬 Homology & Alignment",
        "🌳 Phylogenetics",
        "🧬 Variant Analysis",
        "📑 Data & References",
    ],
)
st.sidebar.divider()
st.sidebar.write("Human MT-CO1")
st.sidebar.write("NCBI Gene ID: 4512")
st.sidebar.write("NC_012920.1:5904..7445")

if page == "🏠 Dashboard":
    st.title("🧬 MitoGene Analyzer")
    st.subheader("In-silico characterization of human mitochondrial MT-CO1")

    c1,c2,c3,c4 = st.columns(4)
    c1.metric("DNA length", f"{len(dna):,} bp")
    c2.metric("Protein length", f"{len(protein):,} aa")
    c3.metric("GC content", f"{dna_stats['GC%']:.2f}%")
    c4.metric("NCBI Gene ID", "4512")

    st.divider()
    st.subheader("Analysis pipeline")
    st.code(
        "NCBI / RefSeq\n"
        "      ↓\n"
        "Reference MT-CO1 sequence\n"
        "      ↓\n"
        "DNA composition → Protein characterization\n"
        "      ↓\n"
        "Homology → Alignment → Phylogeny\n"
        "      ↓\n"
        "Variant consequence analysis",
        language="text",
    )

    st.info(
        "All primary human sequence/protein records are retrieved from NCBI at runtime. "
        "Computed statistics and plots are generated locally by this application."
    )

elif page == "🧬 Gene Information":
    st.title("🧬 Gene Information")
    st.subheader(human["definition"])
    info = pd.DataFrame({
        "Field": ["Gene", "NCBI Gene ID", "Organism", "Reference", "Coordinates", "Protein RefSeq"],
        "Value": [
            "MT-CO1", "4512", "Homo sapiens", "NC_012920.1",
            "5904..7445", "YP_003024028.1"
        ],
    })
    st.dataframe(info, use_container_width=True, hide_index=True)
    st.write(
        "MT-CO1 encodes cytochrome c oxidase subunit I (COX1), a component of "
        "respiratory-chain complex IV in the mitochondrial membrane."
    )
    st.link_button("Open NCBI Gene record", "https://www.ncbi.nlm.nih.gov/gene/4512")

elif page == "🧪 DNA Analysis":
    st.title("🧪 DNA Sequence Analysis")
    c1,c2,c3,c4,c5,c6 = st.columns(6)
    c1.metric("A", dna_stats["A"])
    c2.metric("T", dna_stats["T"])
    c3.metric("G", dna_stats["G"])
    c4.metric("C", dna_stats["C"])
    c5.metric("GC%", f"{dna_stats['GC%']:.2f}")
    c6.metric("AT%", f"{dna_stats['AT%']:.2f}")

    st.subheader("Nucleotide composition")
    comp = pd.DataFrame({"Nucleotide":["A","T","G","C"], "Count":[dna_stats["A"],dna_stats["T"],dna_stats["G"],dna_stats["C"]]})
    st.bar_chart(comp.set_index("Nucleotide"))

    st.subheader("GC sliding-window profile")
    window = st.slider("Window size (bp)", 30, 300, 100, 10)
    step = st.slider("Step size (bp)", 10, 100, 25, 5)
    gc_df = sliding_gc(dna, window, step)
    st.line_chart(gc_df.set_index("Position"))

    with st.expander("Show real reference sequence"):
        st.code(dna, language="text")

    st.download_button(
        "⬇️ Download MT-CO1 FASTA",
        f">MT-CO1_Homo_sapiens|NC_012920.1:5904-7445\n{dna}\n",
        "MT-CO1_Homo_sapiens.fasta",
        "text/plain",
    )

elif page == "🧫 Protein Analysis":
    st.title("🧫 Protein Analysis")
    c1,c2,c3,c4 = st.columns(4)
    c1.metric("Length", f"{len(protein)} aa")
    c2.metric("Molecular weight", f"{prot_stats['Molecular Weight']:.1f} Da")
    c3.metric("pI", f"{prot_stats['Isoelectric Point']:.2f}")
    c4.metric("GRAVY", f"{prot_stats['GRAVY']:.3f}")

    st.write(f"**RefSeq:** {refseq_protein['accession']}")
    st.subheader("Amino-acid composition")
    aa = amino_acid_counts(protein)
    aa_df = pd.DataFrame({"Amino acid": list(aa), "Count": list(aa.values())})
    st.bar_chart(aa_df.set_index("Amino acid"))
    st.dataframe(aa_df, use_container_width=True, hide_index=True)

    st.subheader("Hydrophobicity profile")
    hyd = hydrophobicity(protein)
    st.line_chart(pd.DataFrame({"Residue": range(1,len(hyd)+1), "Kyte-Doolittle": hyd}).set_index("Residue"))

    with st.expander("Show RefSeq protein sequence"):
        st.code(protein, language="text")

    st.download_button(
        "⬇️ Download RefSeq protein FASTA",
        f">YP_003024028.1|COX1_Homo_sapiens\n{protein}\n",
        "YP_003024028.1.fasta",
        "text/plain",
    )

elif page == "🔬 Homology & Alignment":
    st.title("🔬 Homology & Alignment")
    st.write("Homologous mitochondrial CO1/CDS records are retrieved from NCBI and compared locally.")

    if len(homologs) < 2:
        st.warning("Fewer than two homologs were retrieved.")
    else:
        names = list(homologs)
        matrix = pairwise_identity_matrix(homologs)
        st.subheader("Pairwise identity (%)")
        st.dataframe(pd.DataFrame(matrix, index=names, columns=names).round(2), use_container_width=True)

        st.subheader("Select two sequences for global alignment")
        a,b = st.columns(2)
        n1 = a.selectbox("Sequence 1", names)
        n2 = b.selectbox("Sequence 2", names, index=min(1,len(names)-1))
        if n1 != n2:
            result = global_alignment(homologs[n1], homologs[n2])
            st.metric("Identity", f"{result['identity']:.2f}%")
            st.text_area("Alignment", result["alignment"], height=350)
        else:
            st.info("Choose two different sequences.")

        st.subheader("Sequence lengths")
        lengths = pd.DataFrame({"Organism": names, "Length": [len(homologs[n]) for n in names]})
        st.bar_chart(lengths.set_index("Organism"))

elif page == "🌳 Phylogenetics":
    st.title("🌳 Neighbor-Joining Phylogenetics")
    if len(homologs) < 3:
        st.warning("At least three real sequences are needed.")
    else:
        tree, newick = build_nj_tree(homologs)
        fig, ax = plt.subplots(figsize=(10, 7))
        from Bio import Phylo
        Phylo.draw(tree, axes=ax, do_show=False)
        ax.set_title("MT-CO1 Neighbor-Joining Tree")
        st.pyplot(fig)
        plt.close(fig)
        st.download_button("⬇️ Download Newick tree", newick, "MT-CO1_NJ_tree.nwk", "text/plain")
        st.info("This is a distance-based educational tree. It is not presented as a definitive species phylogeny.")

elif page == "🧬 Variant Analysis":
    st.title("🧬 Variant Consequence Analysis")
    st.warning("This module predicts the direct coding consequence of a nucleotide substitution. It does not determine clinical pathogenicity.")
    pos = st.number_input("Reference nucleotide position", 1, len(dna), 1, 1)
    ref = dna[pos-1]
    alt = st.selectbox("Alternative base", ["A","T","G","C"])
    st.write(f"Reference base at position {pos}: **{ref}**")
    if st.button("Analyze substitution"):
        result = analyze_variant(dna, int(pos), alt)
        st.subheader("Result")
        st.json(result)

elif page == "📑 Data & References":
    st.title("📑 Data, Methodology & References")
    st.subheader("Primary data sources")
    st.markdown(
        "- **NCBI Gene:** MT-CO1, Gene ID 4512\n"
        "- **NCBI RefSeq nucleotide:** NC_012920.1, 5904..7445\n"
        "- **NCBI RefSeq protein:** YP_003024028.1\n"
        "- Homologs are retrieved from NCBI mitochondrial records at runtime."
    )
    st.subheader("Methods")
    st.markdown(
        "- Nucleotide counts and GC/AT statistics: calculated locally.\n"
        "- Protein properties: Biopython ProtParam.\n"
        "- Pairwise alignment: Biopython PairwiseAligner.\n"
        "- Phylogeny: Neighbor-Joining from pairwise distances.\n"
        "- Variant consequence: direct sequence/codon comparison using the vertebrate mitochondrial genetic code."
    )
    st.subheader("Important scientific limitation")
    st.info(
        "Database annotations are reported as database records; computed results are generated by this application. "
        "Variant consequence is not a clinical classification."
    )
    st.link_button("NCBI MT-CO1", "https://www.ncbi.nlm.nih.gov/gene/4512")
    st.link_button("NCBI COX1 RefSeq protein", "https://www.ncbi.nlm.nih.gov/protein/YP_003024028.1")
