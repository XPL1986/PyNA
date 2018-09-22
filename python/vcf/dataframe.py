import pandas as pd
from pandas import ExcelWriter
import _pickle as pickle
import os
from python.vcf.vcf import extract, standard


def create_df_list(filenames, location):
    df_list = []
    for filename in filenames:
        data = load_pkl(input_dir(filename, location))
        df_list.append(pd.DataFrame.from_dict(data).set_index(['TAG']))
    return df_list


def concat(df_list):
    """
    Concatenate dataframe for side by side comparision

    :param filenames: [string] of VCF filenames
    :param null: boolean returns only overlapping data
    :return: concatenated dataframe
    """

    blood, tumor, rna = df_list
    unique_fields = ['QUAL', 'MQ', 'DP', 'REF COUNT', 'ALT COUNT', 'VAF']
    new_header = lambda text: {field: field + "_" + text for field in unique_fields}

    shared = ['CHROM', 'POS', 'REF', 'ALT', 'ID', "COMMON", "STANDARD", "GENE", "CDS LENGTH", "FEATURE TYPE", "BIOTYPE",
              "ANNOTATION", "IMPACT"]

    index = pd.concat([blood[shared], tumor[shared], rna[shared]]).drop_duplicates(keep='first')
    new_df = [index, blood.drop(shared + ["EXTRA", "BAIT"], axis=1).rename(columns=new_header("B")),
              tumor.drop(shared + ["EXTRA", "BAIT"], axis=1).rename(columns=new_header("T")),
              rna.drop(shared + ["EXTRA", "BAIT"], axis=1).rename(columns=new_header("R"))]
    df = pd.concat(new_df, axis=1, sort=True)
    _extra_field(df)

    return df


def _extra_field(df, show=True):
    bait = set()
    for genes in df['GENE']:
        for gene in genes.split(";"):
            bait.add(gene)

    score = str(len(bait.intersection(standard))) + "/" + str(len(standard))
    included = ";".join(bait.intersection(standard))
    excluded = ";".join(standard - bait.intersection(standard))
    df["EXTRA"] = [len(df.index.values), score, included, excluded] + [''] * (len(df.index.values) - 4)

    if show:
        print("Score:", score)
        print("Variant Count:", len(df.index.values))
        print("Bait Included:", included)
        print("Bait Excluded:", excluded)
        print("\n")


def df_to_excel(directory, df):
    """Writes dataframe values into CSV files"""

    writer = ExcelWriter(directory)
    df.to_excel(writer, 'Sheet1', index=False)
    writer.save()
    print("Saved to " + directory + "\n")


def input_dir(filename, location):
    if location == 'unfiltered':
        directory = "../../data/" + filename + "/step6/output_ann_dbsnp.vcf"
    elif location == 'gatk':
        directory = "../../data/" + filename + "/snpsift/" + filename + "_annotated_dbsnp_filtered.vcf"
    return directory


def output_dir(filenames, location):

    if location == 'unfiltered':
        path = "/Users/lcwong/Desktop/PyNA/python/output/" + "_".join(filenames) + "/filter_0/"
        if not os.path.exists(path):
            os.makedirs(path)

    elif location == 'gatk':
        path = "/Users/lcwong/Desktop/PyNA/python/output/" + "_".join(filenames) + "/filter_1/"
        if not os.path.exists(path):
            os.makedirs(path)

    return path + "_".join(filenames) + ".xlsx"


def load_pkl(directory):
    try:
        with open(directory[:-4] + ".pkl", 'rb') as f:
            data = pickle.load(f)
    except FileNotFoundError:
        data = extract(directory, save=True)
    return data
