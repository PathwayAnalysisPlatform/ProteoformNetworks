def read_scores(file_name):
    data = pd.read_csv(file_name, sep='\t')
    return data.SCORE
