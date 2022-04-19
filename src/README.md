This is the root folder with all the code for performing Proteoform Network analysis.

The starting point are the jupyter notebook files under the 'python/' directory:

## Set up

- [Get Reactome Graph database as described here.](https://reactome.org/dev/graph-database)

  Suggest using [Neo4j Desktop](https://reactome.org/dev/graph-database/neo4j-desktop) instead of the manual installation.

- [Install Python](https://www.python.org/downloads/)
- Make sure [pip tool](https://pip.pypa.io/en/stable/installation/) for Python is installed:

```
py -m pip install --upgrade pip
```

- Install dependency packages listed in the _requirements.txt_ file. They are necessary for plotting and data handling.

```
py -m pip install -r requirements.txt
```

## Gather reference data:

- Accessioned entities
- Pathways
- Gene - disease associations

## 0. Data overview:

- Number of Genes, Protein and Proteoforms in Reactome.
- Show distribution current knowledge of diversity for each gene (Fig)
- Table with top genes with most proteoforms (Fig)
- Small molecules
- Reactions
- Pathways
- Count how often proteoforms are identical, intermediate and distinct
- Check for protein products of a single gene, what is the distribution of jaccard index of the reactions in which they participate
- Check for proteoform products of a single gene, the distribution of jaccard index values of the reactions in which they participate
- Examples:
  1. Gene with multiple proteoforms

## 1. Network construction:

- Networks by biological process, i.e. per pathway (**Functional**)
- Networks with all known interactions (**Interactome**)
- Networks for a trait of disease (**Module**)

Construct them using a boundary to include only the induced graph from directly related gene products.
Alternatively, construct networks including the first-order neighborhood. Meaning, the induced graph of the closed neighborhood of the directly related nodes.

## 2. Single network characterization

#### Metrics

- 2.1 Size: Number of nodes, links, accessioned entities, small molecules
- 2.2 Node degree: Cases when proteoforms have higher or lower degree with examples
- 2.3 Observability (percolation analysis): Quantify predicting power with threshold proteoform interactome
- 2.4 Connectedness
- 2.5 Articulation points and Bridges
- 2.6 Path length
- 2.7 Average local clustering coefficient

#### Examples

- Cases when proteoform nodes have higher or lower degree with examples. (Fig)
- Examples when proteoforms have distinct set of neighbors, therefore influencing the 'guilty by association' process. (Fig)
- Example of pathways that when modelled with genes vs proteoforms they become disconnected or connected.
- Example of pathways that change much the number of articulation points when changing from genes to proteoforms.
- Examples of pathways that change much the number of bridges when changing from genes to proteoforms.
- Select and plot examples of networks that change their topology when using proteoforms or small molecules
- Metrics of agglomeration of participants of a same disease. The higher the more functional similarity

## 3. Pairwise relationship characterization

#### Metrics

- Network separation index
- Jaccard index
- Overlap index
- Statistics on entities part of the overlap
- Examples of pairs that separate when using proteoforms
- Examples of pairs that share mostly modified proteins
