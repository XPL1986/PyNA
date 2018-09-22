import pandas as pd
import _pickle as pickle
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


def dataframe_to_excel(filename, df):
    output = "/Users/lcwong/Desktop/PyNA/python/output/pipeline_2/" + filename + ".xlsx"
    writer = ExcelWriter(output)
    df.to_excel(writer, 'Sheet1')
    writer.save()


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


def not_common(df, save=False):
    df = df[(df['COMMON'] != 1)]

    print(len(df), len(score(df)), score(df), standard - score(df))
    if save:
        dataframe_to_excel("tumor_freq", df)
    return df


def blood_freq(df, save=False):
    filtered = []
    for index, row in df.iterrows():
        vaf = row['VAF_b']
        dp = row['DP_b']

        if .4 <= vaf <= .6 and dp >= 10:
            filtered.append(index)

    df = df.loc[filtered]
    print(len(df), len(score(df)), score(df), standard - score(df))

    if save:
        dataframe_to_excel("blood_freq", df)
    return df


def tumor_freq(df, save=False):
    filtered = []
    for index, row in df.iterrows():
        vaf = row['VAF_t']
        dp = row['DP_t']

        if pd.isnull(df.loc[index, 'ALT COUNT_t']):
            filtered.append(index)
        elif (vaf <= .1 and dp >= 10) or (vaf >= .9 and dp >= 10):
            filtered.append(index)

    df = df.loc[filtered]
    print(len(df), len(score(df)), score(df), standard - score(df))

    if save:
        dataframe_to_excel("tumor_freq", df)
    return df


def rna_freq(df, save=False):
    filtered = []
    for index, row in df.iterrows():
        vaf = row['VAF_r']
        dp = row['DP_r']

        if pd.isnull(df.loc[index, 'ALT COUNT_r']):
            filtered.append(index)
        elif (vaf <= .1 and dp >= 10) or (vaf >= .9 and dp >= 10):
            filtered.append(index)

    df = df.loc[filtered]
    print(len(df), len(score(df)), score(df), standard - score(df))

    if save:
        dataframe_to_excel("rna_freq", df)
    return df


def run():
    # unfiltered = pd.read_excel(open('/Users/lcwong/Desktop/PyNA/python/output/filter_0/1602B_1602T_CGGA1602.xlsx', 'rb'),
    #                            sheet_name='Sheet1')
    # biotype = protein_coding(unfiltered)
    # dataframe_to_excel("protein_coding", biotype)
    #
    # uncommon = not_common(biotype)
    # dataframe_to_excel("uncommon", uncommon)

    biotype = pd.read_excel(open('/Users/lcwong/Desktop/PyNA/python/output/pipeline_2/protein_coding.xlsx', 'rb'),
                               sheet_name='Sheet1')

    blood = blood_freq(biotype)
    dataframe_to_excel("blood", blood)

    tumor = tumor_freq(blood)
    dataframe_to_excel("tumor_freq", tumor)

    rna = rna_freq(tumor)
    dataframe_to_excel("rna_freq", rna)


run()
