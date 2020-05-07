import os
import sys

LEVELS = ["genes", "proteins", "proteoforms"]
LEVELS_COLOR = {"genes": "teal", "proteins": "teal", "proteoforms": "teal"}

PATH_REACTOME = "resources/reactome/v72/"
FILE_REACTOME_GENES = "genes/all_genes_v72.csv"         # these paths are the suffix of PATH_REACTOME
FILE_REACTOME_PROTEINS = "proteins/all_proteins_v72.csv"
FILE_REACTOME_PROTEOFORMS = "proteoforms/all_proteoforms_v72.csv"
FILE_REACTOME_GENE_INTERACTIONS = "genes/geneInternalEdges.tsv"
FILE_REACTOME_PROTEIN_INTERACTIONS = "proteins/proteinInternalEdges.tsv"
FILE_REACTOME_PROTEOFORM_INTERACTIONS = "proteoforms/proteoformInternalEdges.tsv"
FILE_PROTEOFORMS_SEARCH = "proteoforms/search.tsv"

PATH_RESOURCES = "resources/"
FILE_PATHWAYMATCHER = "PathwayMatcher.jar"
URL_PATHWAYMATCHER = "https://github.com/PathwayAnalysisPlatform/PathwayMatcher/releases/latest/download/PathwayMatcher.jar"

def set_root_wd():
    """Moves to one diretory above the location of the interpreter

    Assumes there is a virtual environment located <repo_root>/venv/
    """
    os.chdir(os.path.dirname(os.path.abspath(sys.executable)) + "\\..\\..")
    print(f"Working directory: {os.getcwd()}")


def get_entity_color(type, level):
    if type == "SimpleEntity":
        return '#8DC7C5'
    else:
        return LEVELS_COLOR[level]