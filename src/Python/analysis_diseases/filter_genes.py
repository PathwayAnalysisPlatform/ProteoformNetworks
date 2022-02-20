import errno
import os
import sys

import pandas as pd

from pathlib import Path

def create_filtered_file(original_file, output_file, significance):
    """

    :param original_file: Original PheGenI file with the gene-trait associations
    :param output_file: Filtered PheGenI file removing the records without genome wide significance
    :param significance: tipicaly
    :return:
    """

    if not Path(original_file).exists():
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), original_file)

    df = pd.read_table(original_file)
    df['P-Value'] = df['P-Value'].astype(float)
    df = df.loc[df['P-Value'] <= 5e-8]
    df.to_csv(output_file, sep="\t", index=False)

    return df

if __name__ == '__main__':
    os.chdir(os.path.dirname(os.path.abspath(sys.executable)) + "\\..\\..\\")
    print(f"Working directory: {os.getcwd()}")
    df = create_filtered_file("resources/PheGenI/PheGenI_Association.txt",
                              "resources/PheGenI/PheGenI_Association_genome_wide_significant.txt",
                              5e-8)

    print(df.shape)
    print(df.columns)
    print(df.dtypes)
