from Bio.Seq import Seq

def analyze_variant(dna, position, alternate):
    if position < 1 or position > len(dna):
        raise ValueError("Position is outside the reference sequence.")
    alternate = alternate.upper()
    if alternate not in {"A","T","G","C"}:
        raise ValueError("Alternative base must be A, T, G or C.")

    ref = dna[position-1]
    if ref == alternate:
        return {
            "position": position,
            "reference": ref,
            "alternate": alternate,
            "consequence": "No nucleotide change",
        }

    mutated = dna[:position-1] + alternate + dna[position:]
    ref_protein = str(Seq(dna).translate(table=2))
    mut_protein = str(Seq(mutated).translate(table=2))

    codon_index = (position-1)//3
    codon_start = codon_index*3
    ref_codon = dna[codon_start:codon_start+3]
    mut_codon = mutated[codon_start:codon_start+3]
    ref_aa = ref_protein[codon_index]
    mut_aa = mut_protein[codon_index]

    if mut_aa == "*":
        consequence = "Stop-gain"
    elif ref_aa == mut_aa:
        consequence = "Synonymous"
    else:
        consequence = "Missense"

    return {
        "position": position,
        "reference": ref,
        "alternate": alternate,
        "reference_codon": ref_codon,
        "alternate_codon": mut_codon,
        "reference_amino_acid": ref_aa,
        "alternate_amino_acid": mut_aa,
        "consequence": consequence,
    }
