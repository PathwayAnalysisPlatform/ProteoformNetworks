import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description='Plot distribution of overlapping scores.')
parser.add_argument('path', metavar="<path>", type=str, help='path to the directory for the plot files')
parser.add_argument('files', metavar='<file_name>', type=str, nargs='+',
                    help='csv file with set overlap scores')

args = parser.parse_args()


def read_scores(file_name):
    data = pd.read_csv(file_name, sep='\t')
    return data.SCORE


# Plot the frequency of each overlap score
# Only considers positive scores, where sets overlap.
def plot(x, score_label='overlap'):
    # Consider only overlapping sets
    x = x[x != 0]
    prob = x.value_counts()
    plt.hist(x, bins=40)
    plt.title(f"Histogram of {score_label} score")
    plt.xlabel('Score')
    plt.ylabel('Frequency')
    plt.show()


# Create the label from the file name.
# All scoring files must follow the naming convention: scores_label.tsv, where 'label' is the name
# of the overlap scoring method
def getScoreLabel(file_name):
    index = file_name.rfind('\\')
    score_label = file_name[(0 if index == -1 else index + 1):]
    score_label = score_label.replace("scores_", "").replace("_", " ").replace(".tsv","")
    return score_label


for file_name in args.files:
    print(f"\n\nPlotting distribution of score file: {file_name} and saving plot to: {args.path}")
    x = read_scores(file_name)
    score_label = getScoreLabel(file_name)
    plot(x, score_label)
