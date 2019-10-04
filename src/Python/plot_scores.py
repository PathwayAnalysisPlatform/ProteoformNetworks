import argparse
import os
import numpy as np
import pandas as pd
from matplotlib.pyplot import plot

import lib.plots
import lib.read_write

parser = argparse.ArgumentParser(description='Plot distribution of overlapping scores.')
parser.add_argument('path', metavar="<path>", type=str, help='path to the directory for the plot files')
parser.add_argument('files', metavar='<file_name>', type=str, nargs='+',
                    help='csv file with set overlap scores')

args = parser.parse_args()


# Main code of the script
for file_name in args.files:
    print(f"\n\nPlotting distribution of score file: {file_name}")
    print(f"Saving plot to: {args.path}")
    x = lib.read_scores(file_name)
    score_label = lib.getScoreLabel(file_name)
    plot(x, score_label)



