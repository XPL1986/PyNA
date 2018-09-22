from python.mutation import translate
from python.mutation.translate import codontable, three_letter


class CodingSequence:
    def __init__(self, transcript_id, strand, nucleotides):
        self.transcript_id = transcript_id
        self.strand = strand
        self.nucleotides = nucleotides
        self.sequence = "".join(list(nucleotides.values()))

    def nucleotide(self, location):
        """Returns nucleotide at given location in index 1"""
        return self.nucleotides[location - 1] if location in self.nucleotides else None

    def mutate(self, location, variant):
        """Mutates the cds at given location in 1 index to given variant"""
        copy = dict()
        for key, value in self.nucleotides.items():
            copy[key] = value
            if key == location - 1:
                copy[key] = variant
        return CodingSequence(self.transcript_id, self.strand, copy)

    def codon(self, location):
        """Finds the corresponding codon at the given location in 1 index"""
        keys = list(self.nucleotides.keys())
        amino_acids = [keys[i:i + 3] for i in range(0, len(keys), 3)]
        for aa in amino_acids:
            if location - 1 in aa:
                codon = "".join([self.nucleotides[i] for i in aa])
        return codon

    def translate(self):
        """Translates sequence into amino acid sequence"""
        return translate.dna_advance(self.sequence)

    @staticmethod
    def single_letter(codon):
        """Translates codon into single letter"""
        return codontable[codon]

    @staticmethod
    def three_letters(codon):
        """Translates codon into three letters"""
        return three_letter[codontable[codon]]


