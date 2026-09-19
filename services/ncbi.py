import time
import requests
from Bio import Entrez, SeqIO
from io import StringIO
from config import NCBI_EMAIL, NCBI_TOOL

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def _params(extra):
    p = {"tool": NCBI_TOOL}
    if NCBI_EMAIL:
        p["email"] = NCBI_EMAIL
    p.update(extra)
    return p

def efetch(db, ident, rettype="fasta", retmode="text", seq_start=None, seq_stop=None):
    params = _params({
        "db": db, "id": ident, "rettype": rettype, "retmode": retmode
    })
    if seq_start is not None:
        params["seq_start"] = seq_start
    if seq_stop is not None:
        params["seq_stop"] = seq_stop
    r = requests.get(f"{BASE}/efetch.fcgi", params=params, timeout=30)
    r.raise_for_status()
    time.sleep(0.35)
    return r.text

def get_human_mtco1():
    fasta = efetch("nuccore", "NC_012920.1", "fasta", "text", 5904, 7445)
    lines = fasta.splitlines()
    seq = "".join(lines[1:]).replace(" ", "").upper()
    return {
        "accession": "NC_012920.1",
        "definition": "MT-CO1 mitochondrially encoded cytochrome c oxidase I [Homo sapiens]",
        "sequence": seq,
    }

def get_protein():
    fasta = efetch("protein", "YP_003024028.1", "fasta", "text")
    record = SeqIO.read(StringIO(fasta), "fasta")
    return {"accession": "YP_003024028.1", "sequence": str(record.seq).upper()}

def _record(accession):
    gb = efetch("nuccore", accession, "gb", "text")
    return SeqIO.read(StringIO(gb), "genbank")

def _find_co1(record):
    for feature in record.features:
        if feature.type != "CDS":
            continue
        gene = " ".join(feature.qualifiers.get("gene", [])).upper()
        product = " ".join(feature.qualifiers.get("product", [])).upper()
        if gene in {"CO1","COX1","COI","MT-CO1","MTCO1"} or "CYTOCHROME C OXIDASE SUBUNIT I" in product or "CYTOCHROME OXIDASE SUBUNIT I" in product:
            return str(feature.extract(record.seq)).upper()
    return None

def get_homologs():
    # NCBI RefSeq complete mitochondrial genomes.
    accessions = {
        "Human": "NC_012920.1",
        "Chimpanzee": "NC_001643.1",
        "Gorilla": "NC_011120.1",
        "Mouse": "NC_005089.1",
        "Rat": "NC_001665.2",
    }
    out = {}
    for name, acc in accessions.items():
        try:
            seq = _find_co1(_record(acc))
            if seq:
                out[name] = seq
        except Exception:
            continue
    return out
