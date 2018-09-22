import re
import _pickle as pickle


def _load_pkl(directory):
    with open(directory, 'rb') as f:
        data = pickle.load(f)
    return data


gene_to_cds_length = _load_pkl("/Users/lcwong/Desktop/PyNA/python/vcf/cds.pkl")

standard = {'PIK3R1', 'PTEN', 'PIK3CG', 'TP53',
            'PTPN11', 'PIK3CA', 'RB1', 'PDGFRA', 'MET',
            'ATRX', 'CDK4', 'EGFR', 'IDH1', 'NF1',
            'CDKN2A', 'MDM4', 'MDM4', 'MDM2', 'CDK6', 'LTBP4'}


def extract(directory, save=False):
    """Extract information from vcf in given directory into dictionary format """

    df_header = ['TAG', 'CHROM', 'POS', 'REF', 'ALT', 'ID', 'COMMON', 'BAIT',
                 'GENE', 'STANDARD', 'CDS LENGTH', 'FEATURE TYPE', 'BIOTYPE',
                 'ANNOTATION', 'IMPACT', 'QUAL', 'MQ', 'DP', 'REF COUNT', 'ALT COUNT', 'VAF']
    df_dict = dict()

    with open(directory) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            else:
                data = _parse_header(line.rstrip().split('\t'))
                for key in df_header:
                    df_dict.setdefault(key, []).append(data[key])

    _extra_field(df_dict)
    if save:
        _save(df_dict, directory)

    return df_dict


def _parse_header(line):
    data = dict()
    vcf_header = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']

    for index, field in enumerate(vcf_header[:6]):
        data[field] = line[index]

    data["TAG"] = ":".join([data['CHROM'], data['POS'], data['REF'], data['ALT']])
    _parse_info(line[7], data)
    _parse_sample(line[9], data)
    return data


def _parse_info(info, data):
    data.setdefault('COMMON', 'N/A')
    for field in info.split(';'):
        if '=' in field:
            key, value = field.split('=')

            if key == 'ANN':
                _parse_ann(value, data)
            elif key == 'MQ':
                data[key] = float(value)
            elif key == 'COMMON':
                data[key] = "FALSE" if "0" in field else "TRUE"

    genes = data["GENE"].split(";")
    bait = set(genes).intersection(standard)
    data["BAIT"] = bait
    data["STANDARD"] = False if len(bait) == 0 else True

    for gene in genes:
        data.setdefault("CDS LENGTH", []).append(str(_calculate_length(gene, gene_to_cds_length)))
    data["CDS LENGTH"] = ";".join(data["CDS LENGTH"])


def _parse_ann(ann, data):
    # ['Allele', 'Annotation', 'Annotation_Impact', 'Gene_Name', 'Gene_ID',
    #  'Feature_Type', 'Feature_ID', 'Transcript_BioType', 'Rank', 'HGVS.c',
    #  'HGVS.p', 'cDNA.pos / cDNA.length', 'CDS.pos / CDS.length', 'AA.pos / AA.length',
    #  'Distance', 'ERRORS / WARNINGS / INFO']
    temp = dict()
    for variant in ann.split(","):
        fields = variant.split("|")
        temp.setdefault("ANNOTATION", set()).add(fields[1])
        temp.setdefault("IMPACT", set()).add(fields[2])
        temp.setdefault("GENE", set()).add(fields[3])
        temp.setdefault("FEATURE TYPE", set()).add(fields[5])
        temp.setdefault("BIOTYPE", set()).add(fields[7])

    for key in ["ANNOTATION", "IMPACT", 'GENE', 'FEATURE TYPE', 'BIOTYPE']:
        data[key] = ";".join(temp[key])


def _parse_sample(sample, data):
    # GT:AD:DP:GQ:PL
    format = re.split("[:,]", sample)
    sample_format = {"REF COUNT": 1, "ALT COUNT": 2, "DP": 3}

    if len(format) == 8:
        for field, index in sample_format.items():
            data[field] = int(format[index])

        vaf = format[2] + "/" + format[3]
        data["VAF"] = round(_convert_to_float(vaf), 3)

    else:
        for field in sample_format:
            data.setdefault(field, 0)
        data.setdefault("VAF", 0)

    return data


def _extra_field(df_dict, show=True):
    bait = set()
    for genes in df_dict['BAIT']:
        bait = bait.union(genes)

    score = str(len(bait)) + "/" + str(len(standard))
    included = ";".join(bait.intersection(standard))
    excluded = ";".join(standard - bait.intersection(standard))
    df_dict["EXTRA"] = [len(df_dict["TAG"]), score, included, excluded] + [''] * (len(df_dict["TAG"]) - 4)

    if show:
        print("Score:", score)
        print("Variant Count:", len(df_dict["TAG"]))
        print("Bait Included:", included)
        print("Bait Excluded:", excluded)


def _calculate_length(gene, gene_to_cds_length):
    """Calculates the average coding sequence length of the gene"""
    try:
        transcripts = gene_to_cds_length[gene]
    except KeyError:
        transcripts = []

    lengths = []
    for transcript in transcripts:
        lengths.append(transcripts[transcript]["length"])

    length = round(sum(lengths) / float(len(lengths)) if len(lengths) != 0 else 0)
    return length


def _convert_to_float(frac_str):
    """Converts a string into a float"""

    try:
        return float(frac_str)
    except ValueError:
        num, denom = frac_str.split('/')
        try:
            leading, num = num.split(' ')
            whole = float(leading)
        except ValueError:
            whole = 0

        if float(denom) == 0:
            return 0

        fraction = float(num) / float(denom)
        return whole - fraction if whole < 0 else whole + fraction


def _save(data, directory):
    """Saves the dictionary pickle files"""
    directory = directory[:-4] + ".pkl"
    with open(directory, 'wb') as f:
        pickle.dump(data, f, -1)
    print("Saved to " + directory + "\n")
