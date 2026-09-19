# MitoGene Analyzer — Advanced

An educational Streamlit application for in-silico characterization of human mitochondrial MT-CO1 using live NCBI/RefSeq records.

## Real data

Human reference:
- NCBI Gene: 4512
- Nucleotide: NC_012920.1, 5904..7445
- Protein: YP_003024028.1

The application retrieves these records from NCBI E-utilities at runtime.

## Run

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
streamlit run app.py
```

## NCBI request identification

Before public deployment, set `NCBI_EMAIL` in `config.py` to the email of the software developer. Keep request volume low and respect NCBI E-utilities policies.

## Scientific scope

The app calculates sequence statistics, protein properties, pairwise sequence identity, a Neighbor-Joining tree, and direct coding consequences of simulated substitutions.

Variant consequence is not a clinical pathogenicity classification.
