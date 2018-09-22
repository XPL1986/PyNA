import os
import pysam
from pandas import ExcelWriter
import pandas as pd

references = {"GRCh37": [], "NCBI36": [], "other": [], "error": []}

def sort_bam(path):
    for filename in os.listdir(path):
        if filename.endswith(".bam"):
            print(filename)
            try:
                bam = pysam.AlignmentFile(path + filename, "rb")
                chr1_length = bam.header["SQ"][0]['LN']

                if chr1_length == 249250621:
                    references["GRCh37"].append(filename)
                elif chr1_length == 247249719:
                    references["NCBI36"].append(filename)
                else:
                    references["other"].append(filename)
            except OSError:
                references["error"].append(filename)


sort_bam("/Volumes/qmu/NG2016/rGBM_DNA/")
sort_bam("/Volumes/qmu/NG2016/rGBM_RNA/")

base = len(references[max(references, key=lambda a: len(references[a]))])
for key, value in references.items():
    references[key] = references[key] + ([""] * (base - len(references[key])))

writer = ExcelWriter("sorted.xlsx")
df = pd.DataFrame.from_dict(references)
df.to_excel(writer, 'Sheet1', index=False)
writer.save()

