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
    index = file_name.rfind('/')
    score_label = file_name[(0 if index == -1 else index + 1):]
    score_label = score_label.replace("scores_", "").replace("_", " ").replace(".tsv", "")
    return score_label