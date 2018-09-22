import re

from python.gtf import request
from python.mutation import fasta, rna
from python.mutation.translate import codontable, three_letter


def _find_cds(location):
    """Finds associated transcript coding sequence at given chr:coord"""

    location_to_gene = request.load('../gtf/dictionaries/location.json')
    features = request.location_lookup(location, location_to_gene)
    cds_list = [feature["transcript_id"] for feature in features if feature["feature"] == "CDS"]
    return cds_list


def _complement(sequence):
    """returns the complementary sequence of given sequence"""
    base_complement = {"A": "T", "T": "A", "G": "C", "C": "G"}

    reverse = []
    for nucleotide in sequence:
        if nucleotide in base_complement:
            reverse.append(base_complement[nucleotide])
        else:
            print("contains N or intronic regions.")
            reverse.append(nucleotide)

    reverse = "".join(reverse)[::-1]
    return reverse


def _create_cds(cds_list, fa):
    if len(cds_list) == 0:
        return "Not located in any known coding sequence"

    transcript_to_cds = request.load('../gtf/dictionaries/cds.json')
    cds = []

    for transcript_id in cds_list:
        sequence, location = [], []
        strand = None
        for exon in transcript_to_cds[transcript_id]["cds"]:
            exon = re.split("[:[,)/]", exon)
            start, end, strand = int(exon[2]), int(exon[3]), exon[5]

            if strand == "+":
                sequence.append(fasta.nucleotides(fa, start, end))
                location += [i for i in range(start, end)]
            elif strand == "-":
                sequence.append(_complement(fasta.nucleotides(fa, start, end)))
                location += [i for i in range(end - 1, start - 1, -1)]

        sequence = "".join(sequence)
        nucleotides = {location[i]: sequence[i] for i in range(len(location))}
        cds.append(rna.CodingSequence(transcript_id, strand, nucleotides))

    return cds


def snp(seqname, location, variant):
    coordinate = seqname + ":" + str(location)
    fa = fasta.load(seqname + ".fa")
    cds_list = _create_cds(_find_cds(coordinate), fa)

    for cds in cds_list:
        print("Transcript ID:", cds.transcript_id)
        print("Original Sequence:", cds.sequence)
        print("Original Protein:", cds.translate())
        print("Nucleotide Counts:", len(cds.sequence))
        print("Amino Acid Counts:", len(cds.translate()) - 1)

        base_complement = {"A": "T", "T": "A", "G": "C", "C": "G"}

        mutated = cds.mutate(location, variant) if cds.strand == "+" else cds.mutate(location, base_complement[variant])
        ref = cds.nucleotide(location) if cds.strand == "+" else base_complement[cds.nucleotide(location)]

        print(seqname + ":g." + str(location) + ref + ">" + variant)
        print(codontable[cds.codon(location)] + "(" + three_letter[codontable[cds.codon(location)]] + ")" + ">"
              + codontable[mutated.codon(location)] + "(" + three_letter[codontable[mutated.codon(location)]] + ")")

        print("Original Sequence:", mutated.sequence)
        print("Original Protein:", mutated.translate())
        print("Nucleotide Counts:", len(mutated.sequence))
        print("Amino Acid Counts:", len(mutated.translate()) - 1)
        print("Genome Build: GRCh38 / hg38")


# snp("chr10", 87933148, "A")
snp("chr2", 208248418, "T")
