import pandas as pd
from pandas import ExcelWriter

data = dict()

standard = {'PIK3R1', 'PTEN', 'PIK3CG', 'TP53',
            'PTPN11', 'PIK3CA', 'RB1', 'PDGFRA', 'MET',
            'ATRX', 'CDK4', 'EGFR', 'IDH1', 'NF1',
            'CDKN2A', 'MDM4', 'MDM4', 'MDM2', 'CDK6', 'LTBP4'}


def score(df):
    bait = set()
    for genes in df["GENE"]:
        for gene in str(genes).split(";"):
            bait.add(gene)

    intersection = standard.intersection(bait)
    return intersection


def tumor_freq(df, t, dp):
    df = df[(df['VAF_t'] >= t) & (df['DP_t'] >= dp)]
    return df


def rna_freq(df, r, dp):
    df = df[(df['VAF_r'] >= r) & (df['DP_r'] >= dp)]
    return df


def length(df, l):
    df = df[df['LENGTH'].apply(lambda x: any(int(i) < l for i in x.split(";")))]
    return df


def optimize(df):
    tumor_range = [i for i in range(1, 11, 1)]
    rna_range = [i for i in range(1, 92, 1)]
    length_range = [i for i in range(6500, 2000, -300)]

    # tumor filter
    for t in tumor_range:
        df_t = tumor_freq(df, t/100, 25)
        if not truth(df_t):
            break

        # rna filter
        for r in rna_range:
            df_r = rna_freq(df_t, r/100, 25)

            if not truth(df_r):

                break

            # length filter
            for l in length_range:
                df_l = length(df_r, l)

                if not truth(df_l):
                    break

                write(df_l, t, r, l)

    df = pd.DataFrame.from_dict(data)
    return df


def write(df, t, r, l):
    data.setdefault("tumor_freq", []).append(t)
    data.setdefault("rna_freq", []).append(r)
    data.setdefault("length", []).append(l)
    data.setdefault("variants", []).append(len(df))

    bait = score(df)
    data.setdefault("score", []).append(len(bait))
    data.setdefault("bait", []).append(bait)
    return df


def truth(df):
    return len(score(df)) > 0

def dataframe_to_excel(filename, df):
    """
    Writes dataframe values into CSV files

    :param filename:
    :param df: dataframe
    :return: None
    """
    output = filename + ".xlsx"
    writer = ExcelWriter(output)
    df.to_excel(writer, 'Sheet1')
    writer.save()

filtered = pd.read_excel(open('/Users/lcwong/Desktop/PyNA/python/output/pipeline_1/annotation.xlsx', 'rb'), sheet_name='Sheet1')
dataframe_to_excel("tune", optimize(filtered))
