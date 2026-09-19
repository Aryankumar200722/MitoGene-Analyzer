from Bio.SeqUtils.ProtParam import ProteinAnalysis

KD = {
    "I":4.5,"V":4.2,"L":3.8,"F":2.8,"C":2.5,"M":1.9,"A":1.8,"G":-0.4,
    "T":-0.7,"S":-0.8,"W":-0.9,"Y":-1.3,"P":-1.6,"H":-3.2,"E":-3.5,
    "Q":-3.5,"D":-3.5,"N":-3.5,"K":-3.9,"R":-4.5
}

def protein_stats(seq):
    p = ProteinAnalysis(seq)
    return {
        "Molecular Weight": p.molecular_weight(),
        "Aromaticity": p.aromaticity(),
        "Instability Index": p.instability_index(),
        "Isoelectric Point": p.isoelectric_point(),
        "GRAVY": p.gravy(),
    }

def amino_acid_counts(seq):
    return {aa: seq.count(aa) for aa in "ACDEFGHIKLMNPQRSTVWY"}

def hydrophobicity(seq, window=19):
    vals = []
    half = window // 2
    for i in range(len(seq)):
        left = max(0, i-half)
        right = min(len(seq), i+half+1)
        chunk = seq[left:right]
        vals.append(sum(KD.get(x,0) for x in chunk)/len(chunk))
    return vals
