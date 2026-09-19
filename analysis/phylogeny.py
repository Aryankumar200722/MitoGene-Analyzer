from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo
from analysis.alignment import pairwise_identity

def build_nj_tree(sequences):
    names = list(sequences)
    lower = []
    for i in range(len(names)):
        row = []
        for j in range(i+1):
            if i == j:
                row.append(0)
            else:
                identity = pairwise_identity(sequences[names[i]], sequences[names[j]])
                row.append(1 - identity/100)
        lower.append(row)
    dm = DistanceMatrix(names, lower)
    tree = DistanceTreeConstructor().nj(dm)
    import io
    buf = io.StringIO()
    Phylo.write(tree, buf, "newick")
    return tree, buf.getvalue()
