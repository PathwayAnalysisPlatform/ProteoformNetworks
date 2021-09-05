import os
import sys

import networkx as nx
from bokeh.io import show

from config import proteoforms, with_unique_sm
from lib.networks import create_pathway_interaction_network
from visualization.visualize_single_network import plot_interaction_network, Coloring


def main():
    pathway1 = "R-HSA-9613829"  # Regulation of glycolysis by fructose 2,6-bisphosphate metabolism
    os.chdir(os.path.dirname(os.path.abspath(sys.executable)) + "\\..\\..")
    print(f"Working directory: {os.getcwd()}")

    g = create_pathway_interaction_network(pathway1, proteoforms, with_unique_sm, "resources/pathway_networks/")

    # Calculate articulation points
    list()


    # Calculate bridges

    # Mark edges

    # Update the graph file

    p = plot_interaction_network(g, coloring=Coloring.ENTITY_TYPE, plot_width=600, plot_height=500,
                                title="Test Plot Bridges and Articulation Points",
                                legend_location="right")
    show(p)
    print("Finished")


if __name__ == '__main__':
    main()