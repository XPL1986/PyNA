from itertools import takewhile

codontable = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}

three_letter = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp",
    "C": "Cys", "E": "Glu", "Q": "Gln", "G": "Gly",
    "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
    "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser",
    "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "_": "Stop"
}


def dna_simple(sequence):
    """
    :param sequence: sequence of nucleotide
    :return: simple translation
    """
    codons = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]
    aminoacids = [codontable[codon] for codon in codons if len(codon) == 3]
    return "-".join(aminoacids)


# Source: https://stackoverflow.com/questions/19521905/translation-dna-to-protein

def dna_advance(sequence, codontable=codontable, stop_codons=('TAA', 'TGA', 'TAG')):
    start = sequence.find('ATG')

    # Take sequence from the first start codon
    trimmed_sequence = sequence[start:]

    # Split it into triplets
    codons = [trimmed_sequence[i:i + 3] for i in range(0, len(trimmed_sequence), 3)]

    # Take all codons until first stop codon
    coding_sequence = takewhile(lambda x: x not in stop_codons and len(x) == 3, codons)

    # Translate and join into string
    protein_sequence = ''.join([codontable[codon] for codon in coding_sequence])
    # This line assumes there is always stop codon in the sequence
    return "{0}_".format(protein_sequence)
