import pandas as pd
import os

scores_genes = pd.read_csv("../../reports/scores_genes_jaccard_similarity.tsv", sep='\t')
scores_proteins = pd.read_csv("../../reports/scores_proteins_jaccard_similarity.tsv", sep='\t')
scores_proteoforms = pd.read_csv("../../reports/scores_proteoforms_jaccard_similarity.tsv", sep='\t')

scores = pd.DataFrame(scores_genes, scores_proteins, scores_proteoforms)

print(scores.info())
