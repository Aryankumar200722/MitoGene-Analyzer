from Bio.Align import PairwiseAligner

def global_alignment(a, b):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aln = aligner.align(a,b)[0]
    s1, s2 = str(aln[0]), str(aln[1])
    compared = 0
    matches = 0
    for x,y in zip(s1,s2):
        if x != "-" and y != "-":
            compared += 1
            matches += (x == y)
    identity = matches/compared*100 if compared else 0
    return {"identity": identity, "alignment": str(aln)}

def pairwise_identity(a,b):
    return global_alignment(a,b)["identity"]

def pairwise_identity_matrix(sequences):
    names = list(sequences)
    return {
        n1: {n2: pairwise_identity(sequences[n1], sequences[n2]) for n2 in names}
        for n1 in names
    }
