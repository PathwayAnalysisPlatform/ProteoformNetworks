import os
import sys

from bokeh.palettes import Colorblind

LEVELS = ["genes", "proteins", "proteoforms"]

# LEVELS_COLOR = {"genes": "#67A9CF", "proteins": "#F1A340", "proteoforms": "#7FBF7B"}
LEVELS_COLOR = {"genes": Colorblind[4][0], "proteins": Colorblind[4][1], "proteoforms": Colorblind[4][3]}
PROTEOFORMS_SUBGRAPHS_COLOR = {"mm": Colorblind[6][5], "uu": Colorblind[6][2], "um": Colorblind[6][1]}

PATH_REACTOME = "resources/reactome/"
FILE_REACTOME_GENES = "genes_vertices.tsv"  # these paths are the suffix of PATH_REACTOME
FILE_REACTOME_PROTEINS = "proteins_vertices.tsv"
FILE_REACTOME_PROTEOFORMS = "proteoforms_vertices.tsv"
FILE_REACTOME_GENE_INTERACTIONS = "genes_interactions.tsv"
FILE_REACTOME_PROTEIN_INTERACTIONS = "proteins_interactions.tsv"
FILE_REACTOME_PROTEOFORM_INTERACTIONS = "proteoforms_interactions.tsv"
FILE_PROTEOFORMS_SEARCH = "proteoforms/search.tsv"

PATH_RESOURCES = "resources/"
FILE_PATHWAYMATCHER = "PathwayMatcher.jar"
URL_PATHWAYMATCHER = "https://github.com/PathwayAnalysisPlatform/PathwayMatcher/releases/latest/download/PathwayMatcher.jar"

GRAPHS_PATH = "resources/Reactome/"
MAPPING_FILE = "mapping_proteins_to_level.tsv"


def set_root_wd():
    """Moves to one diretory above the location of the interpreter

    Assumes there is a virtual environment located <repo_root>/venv/
    """
    os.chdir(os.path.dirname(os.path.abspath(sys.executable)) + "\\..\\..")
    print(f"Working directory: {os.getcwd()}")


def get_entity_color(type, level):
    if type == "SimpleEntity":
        return '#BDBDBD'
    else:
        return LEVELS_COLOR[level]


COLOR_IO = "black"
COLOR_CO = "orange"
COLOR_RO = "red"
COLOR_CC = "dodgerblue"

EDGES_WEIGHT_IN_REACTION = 1
EDGES_WEIGHT_IN_COMPLEX = 50
