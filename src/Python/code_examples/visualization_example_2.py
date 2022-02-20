import pandas as pd
from bokeh.io import show, curdoc
from bokeh.layouts import layout
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import DataTable, TableColumn, Div

from config import LEVELS
from visualization.visualize_two_networks import create_overlap_graph, plot_overlap_graph

trait1 = "Adiponectin"
trait2 = "Thyrotropin"
path_to_modules = "../../reports/All_modules/"
path_to_figures = "../../figures/overlap_analysis/"

overlap_data = pd.read_table("../../reports/overlap_similarity_score_variation_examples.tsv")
overlap_data.head()

data = dict(
    trait1=overlap_data["TRAIT1"],
    trait2=overlap_data["TRAIT2"],
    score=overlap_data["PROTEINS_TO_PROTEOFORMS"]
)
source = ColumnDataSource(data)

columns = [
    TableColumn(field="trait1", title="Trait1"),
    TableColumn(field="trait2", title="Trait2"),
    TableColumn(field="score", title="Score"),
]
data_table = DataTable(source=source, columns=columns, width=600, height=280)
show(data_table)


def callback(attrname, old, new):
    print("Attribute changed: ", attrname)
    print("Old: ", old)
    print("New: ", new)
    selectionIndex = source.selected.indices[0]
    print("you have selected the row numr " + str(selectionIndex) + " which is "
          + overlap_data["TRAIT1"].iloc[selectionIndex] + " -- " + overlap_data["TRAIT2"].iloc[selectionIndex])
    trait1 = overlap_data["TRAIT1"].iloc[selectionIndex]
    trait2 = overlap_data["TRAIT2"].iloc[selectionIndex]


source.selected.on_change("indices", callback);

title = f"<p style=\"font-weight:bold;text-align:center;font-size:22px;width:1800px;\">"             f"<span style=\"color:green;\">{trait1}</span> with "             f"<span style=\"color:blue;\">{trait2}</span>"             f"</p>"

graphs_complete = {level: create_overlap_graph(trait1, trait2, level, path_to_modules) for level in LEVELS}
graphs_interface = {level: create_overlap_graph(trait1, trait2, level, path_to_modules, only_interface=True)
                    for level in LEVELS}

figures_complete_modules = [plot_overlap_graph(level, graph) for level, graph in graphs_complete.items()]
figures_interfaces = [plot_overlap_graph(level, graph) for level, graph in graphs_interface.items()]

l = layout(
    [[data_table],
     [Div(text=f"{title}")],
     figures_complete_modules,
     figures_interfaces,
     ])

show(l)
curdoc().add_root(l)
curdoc().title = "Module pairs visualization"
