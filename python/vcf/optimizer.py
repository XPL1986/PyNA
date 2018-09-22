import pandas as pd
import re
import json
from pandas import ExcelWriter

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


def blood_freq(df, save=False):
    filtered = []
    for index, row in df.iterrows():
        ad = row['ALT COUNT_b']
        dp = row['DP_b']

        if pd.isnull(df.loc[index, 'ALT COUNT_b']):
            filtered.append(index)
        elif (ad == 0) or (ad == 1 and dp > 25) or (ad == 2 and dp > 50):
            filtered.append(index)

    df = df.loc[filtered]
    print(len(df), len(score(df)), score(df), standard - score(df))

    if save:
        dataframe_to_excel("blood_freq", df)
    return df




def protein_coding(df, save=False):
    filtered = []
    for index, row in df.iterrows():

        if not pd.isnull(df.loc[index, 'BIOTYPE']):
            if "protein_coding" in row['BIOTYPE'].split(";"):
                filtered.append(index)

    df = df.loc[filtered]
    print(len(df), len(score(df)), score(df), standard - score(df))
    if save:
        dataframe_to_excel("biotype", df)
    return df


def impact(df, save=False):
    df = df[df['IMPACT'].apply(lambda x: False if len(x.split(";")) == 1 and x.split(";")[0] == 'LOW' else True)]
    print(len(df), len(score(df)), score(df), standard - score(df))

    if save:
        dataframe_to_excel("impact", df)
    return df


def annotation(df, save=False):
    parameters = pd.read_excel(open('/Users/lcwong/Desktop/PyNA/python/vcf/annotations/annotations.xlsx', 'rb'),
                               sheet_name='Sheet1')

    important = set(parameters['Important'].tolist())
    high = set(parameters['High Only'].tolist())
    is_high = lambda line: "HIGH" in line['IMPACT'].split(";")

    filtered = []
    for index, row in df.iterrows():
        ann = set(re.split('[&;]+', row['ANNOTATION']))
        if len(ann.intersection(important)) != 0:
            filtered.append(index)
        elif len(ann.intersection(high)) != 0 and is_high(row):
            filtered.append(index)

    df = df.loc[filtered]

    print(len(df), len(score(df)), score(df), standard - score(df))
    if save:
        dataframe_to_excel("annotation", df)
    return df


def tumor_freq(df, t, dp, save=False):
    df = df[(df['VAF_t'] >= t) & (df['DP_t'] >= dp)]

    print(len(df), len(score(df)), score(df), standard - score(df))

    if save:
        dataframe_to_excel("tumor_freq", df)
    return df



def rna_freq(df, r, dp, save=False):
    df = df[(df['VAF_r'] >= r) & (df['DP_r'] >= dp)]

    print(len(df), len(score(df)), score(df), standard - score(df))

    if save:
        dataframe_to_excel("rna_freq", df)
    return df

# def rna_freq(df, r, dp, save=False):
#     df = df[(df['VAF_r'] <= .1) & (df['DP_r'] >= 10) | (df['VAF_r'] >= .9) & (df['DP_r'] >= 10)]
#
#     print(len(df), len(score(df)), score(df), standard - score(df))
#
#     if save:
#         dataframe_to_excel("rna_freq", df)
#     return df


def length(df, l, save=False):
    df = df[df['LENGTH'].apply(lambda x: any(int(i) < l for i in x.split(";")))]

    print(len(df), len(score(df)), score(df), standard - score(df))
    if save:
        dataframe_to_excel("length", df)
    return df





def dataframe_to_excel(filename, df):
    """
    Writes dataframe values into CSV files

    :param filename:
    :param df: dataframe
    :return: None
    """
    output = "output/backup/" + filename + ".xlsx"
    writer = ExcelWriter(output)
    df.to_excel(writer, 'Sheet1')
    writer.save()


def tune(t, r, l):
    with open("vcf/pipeline_1/gatk.json") as f:
        df = pd.DataFrame.from_dict(json.load(f))
        print(len(df), len(score(df)), score(df), standard - score(df))
        f.close()

    df = impact(df, save=False)
    df = annotation(df, save=False)

    df = tumor_freq(df, t, save=False)
    df = rna_freq(df, r, save=False)
    df = length(df, l, save=False)

    dataframe_to_excel("trial_" + "_".join([str(t), str(r), str(l)]), df)
    print(df)


if __name__ == '__main__':
    # df = pd.read_excel(open('/Users/lcwong/Desktop/PyNA/python/output/pipeline_1/Unfiltered.xlsx', 'rb'), sheet_name='Sheet1')
    # dataframe_to_excel("/Users/lcwong/Desktop/PyNA/python/output/pipeline_1/rna_freq.xlsx", blood_freq(df, save=False))
    tune(1, 1, 1)
