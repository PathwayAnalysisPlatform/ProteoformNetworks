# Map gene ids
import pandas as pd

from lib.conversions import create_gene_to_protein_mapping

path_file_genes = "../../../resources/Reactome/v70/Genes/all_genes_v70.csv"

create_gene_to_protein_mapping(path_file_genes, "../../../resources/UniProt/", "mapping_genes_proteins.csv", 1000)