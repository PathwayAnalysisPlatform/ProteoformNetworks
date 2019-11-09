import os
import re

import networkx as nx
import pandas as pd

from lib.download import download_if_not_exists

LEVELS = {"genes", "proteins", "proteoforms"}


def read_scores(file_name):
    data = pd.read_csv(file_name, sep='\t')
    return data.SCORE


def get_trait_file_name(trait):
    return re.sub("[ ,]", "_", trait.replace("\"", ""))


def get_graph(trait, level, path_to_modules):
    """ Create a networkx graph instance for a trait module at the level specified """
    if level not in LEVELS:
        raise ValueError("level must be one of %r." % LEVELS)

    G = nx.Graph()

    # Traverse file with module members. Get set of members for the trait
    print("Reading members of module ", trait)
    prefix = path_to_modules + get_trait_file_name(trait) + "_" + level
    path_file_module_vertices = prefix + "_vertices.tsv"

    with open(path_file_module_vertices) as file_vertices:
        file_vertices.readline()  # Read header
        line = file_vertices.readline()
        while line:
            G.add_node(line.strip())
            line = file_vertices.readline()

    # Traverse file with interactions of the level. Get the set of edges for this trait
    print("Reading interactions in module ", trait)
    path_file_module_edges = prefix + "_edges.tsv"

    with open(path_file_module_edges) as file_edges:
        file_edges.readline()  # Read header
        line = file_edges.readline()
        while line:
            columns = line.split("\t")
            G.add_edge(columns[0], columns[1].strip())
            line = file_edges.readline()

    return G


def create_pathwaymatcher_files(path_reactome,
                                file_reactome_genes, file_reactome_proteins, file_reactome_proteoforms,
                                file_reactome_gene_interactions,
                                file_reactome_protein_interactions,
                                file_reactome_proteoform_interactions,
                                path_pathwaymatcher, file_pathwaymatcher, url_pathwaymatcher):
    """Create protein interaction network of all Reactome by mapping with PathwayMatcher.
    Also, creates the search files for genes, proteins and proteoforms, to be used as the sets to create the modules."""

    assert os.path.exists(path_reactome + file_reactome_genes), f"Missing file: {path_reactome}{file_reactome_genes}"
    assert os.path.exists(
        path_reactome + file_reactome_proteins), f"Missing file: {path_reactome}{file_reactome_proteins}"
    assert os.path.exists(
        path_reactome + file_reactome_proteoforms), f"Missing file: {path_reactome}{file_reactome_proteoforms}"

    if not os.path.exists(path_reactome + file_reactome_gene_interactions) \
            or not os.path.exists(path_reactome + file_reactome_protein_interactions) \
            or not os.path.exists(path_reactome + file_reactome_proteoform_interactions):

        # Download PathwayMatcher executable
        download_if_not_exists(path_pathwaymatcher, file_pathwaymatcher, url_pathwaymatcher, 'PathwayMatcher')

        print("Creating gene interaction file...")
        if not os.path.exists(path_reactome + file_reactome_gene_interactions):
            command = f"java -jar {path_pathwaymatcher}{file_pathwaymatcher} match-genes -i {path_reactome}{file_reactome_genes} -o {path_reactome}Genes/ -g"
            os.system(command)

        print("Creating protein interaction file...")
        if not os.path.exists(path_reactome + file_reactome_protein_interactions):
            command = f"java -jar {path_pathwaymatcher}{file_pathwaymatcher} match-uniprot -i {path_reactome}{file_reactome_proteins} -o {path_reactome}Proteins/ -g"
            os.system(command)

        print("Creating proteoform interaction file...")
        if not os.path.exists(path_reactome + file_reactome_proteoform_interactions):
            command = f"java -jar {path_pathwaymatcher}{file_pathwaymatcher} match-proteoforms -i {path_reactome}{file_reactome_proteoforms} -o {path_reactome}Proteoforms/ -g"
            os.system(command)

        # Remove extra PathwayMatcher result files
        extra_files = ["Genes/analysis.tsv", "Proteins/analysis.tsv", "Proteoforms/analysis.tsv",
                       "Genes/geneExternalEdges.tsv", "Proteins/proteinExternalEdges.tsv",
                       "Proteoforms/proteoformExternalEdges.tsv",
                       "Genes/geneVertices.tsv", "Proteins/proteinVertices.tsv", "Proteoforms/proteoformVertices.tsv",
                       "Genes/search.tsv", "Proteins/search.tsv", "Proteoforms/search.tsv"]
        for file in extra_files:
            if os.path.exists(f"{path_reactome}{file}"):
                os.remove(f"{path_reactome}{file}")

    print(f"PathwayMatcher files READY")
