import pandas as pd
import re
import os
import argparse
from pandas import ExcelWriter
from python.vcf.dataframe import create_df_list, concat, df_to_excel, output_dir
from python.vcf.vcf import standard


def protein_coding(df, save, bundle):
    passed = []
    for index, row in df.iterrows():
        if not pd.isnull(df.loc[index, 'BIOTYPE']):
            if "protein_coding" in row['BIOTYPE'].split(";"):
                passed.append(index)
    df = df.loc[passed]

    if save:
        dataframe_to_excel("1_Biotype", df, bundle)

    return df


def uncommon(df, save, bundle):
    df = df[(df['COMMON'] != 'TRUE')]

    if save:
        dataframe_to_excel("2_Uncommon", df, bundle)
    return df


def blood_freq(df, save, bundle):
    filtered = []
    for index, row in df.iterrows():
        ad, dp = row['ALT COUNT_B'], row['DP_B']

        if pd.isnull(df.loc[index, 'ALT COUNT_B']):
            filtered.append(index)
        elif (ad == 0) or (ad == 1 and dp > 25) or (ad == 2 and dp > 50):
            filtered.append(index)
    df = df.loc[filtered]

    if save:
        dataframe_to_excel("3_Blood_freq", df, bundle)
    return df


def gatk_filter(df, save, bundle):
    df = blood_freq(uncommon(protein_coding(df, False, bundle), False, bundle), False, bundle)

    if save:
        dataframe_to_excel("4_GATK", df, bundle)
    return df


def impact(df, save, bundle):
    df = df[df['IMPACT'].apply(lambda x: False if len(x.split(";")) == 1 and x.split(";")[0] == 'LOW' else True)]

    if save:
        dataframe_to_excel("5_Impact", df, bundle)
    return df


def annotation(df, save, bundle):
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

    if save:
        dataframe_to_excel("6_Annotation", df, bundle)
    return df


def tumor_freq(df, t, dp, save, bundle):
    df = df[(df['VAF_T'] >= t) & (df['DP_T'] >= dp)]

    if save:
        dataframe_to_excel("7_Tumor_freq", df, bundle)
    return df


def rna_freq(df, r, dp, save, bundle):
    df = df[(df['VAF_R'] >= r) & (df['DP_R'] >= dp)]

    if save:
        dataframe_to_excel("8_RNA_freq", df, bundle)
    return df


def length(df, l, save, bundle):
    df = df[df['CDS LENGTH'].apply(lambda x: any(int(i) <= l for i in x.split(";")))]

    if save:
        dataframe_to_excel("9_Length", df, bundle)
    return df


def info(df):
    bait = set()
    for genes in df["GENE"]:
        for gene in str(genes).split(";"):
            bait.add(gene)

    score = str(len(bait.intersection(standard))) + "/" + str(len(standard))
    included = ";".join(bait.intersection(standard))
    excluded = ";".join(standard - bait.intersection(standard))

    print("Score:", score)
    print("Variant Count:", len(df.index.values))
    print("Bait Included:", included)
    print("Bait Excluded:", excluded)
    print("\n")

    return df


def dataframe_to_excel(filename, df, bundle):
    path = "/Users/lcwong/Desktop/PyNA/python/output/" + "_".join(bundle) + "/pipeline_1/"
    if not os.path.exists(path):
        os.makedirs(path)

    output = path + filename + ".xlsx"
    writer = ExcelWriter(output)
    df.to_excel(writer, 'Sheet1', index=False)
    writer.save()
    print("Saved to " + output)


def run(bundle):
    # filter_0 = concat(create_df_list(bundle, 'unfiltered'))
    # df_to_excel(output_dir(bundle, 'unfiltered'), filter_0)
    #
    # filter_1 = info(protein_coding(filter_0, True, bundle))
    # filter_2 = info(uncommon(filter_1, True, bundle))
    # filter_3 = info(blood_freq(filter_2, True, bundle))

    gatk = concat(create_df_list(bundle, 'gatk'))
    # df_to_excel(output_dir(bundle, 'gatk'), gatk)

    filter_4 = info(gatk_filter(gatk, False, bundle))
    filter_5 = info(impact(filter_4, False, bundle))
    filter_6 = info(annotation(filter_5, True, bundle))
    filter_7 = info(tumor_freq(filter_6, .05, 25, True, bundle))
    filter_8 = info(rna_freq(filter_7, .3, 25, True, bundle))
    filter_9 = info(length(filter_8, 3000, True, bundle))

run(["R056.DNA.0", "R056.DNA.2", "R056.RNA.2"])
# run(["R060.DNA.0", "R060.DNA.2", "R060.RNA.2"])
