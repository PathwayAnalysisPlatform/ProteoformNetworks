This is the root folder with all the code for performing Proteoform Network analysis.

The starting point are the jupyter notebook files under the 'python/' directory:

### 0. Pre-analysis

#### 0.0 Set up

- Install Neo4j
- Run python

#### 0.1 Obtain all reference data:

- Accessioned entities
- Pathways
- Gene - disease associations

#### 0.2 Reference data overview:

- Genes
- Proteins
- Proteoforms
- Show distribution current knowledge of diversity for each gene (Fig)
- Table with top genes with most proteoforms (Fig)
- Small molecules
- Reactions
- Pathways

### 1. Network construction:

1.1 Interactomes
1.2 Pathways into gene, protein and proteoform interaction networks
1.3 Disease networks with gene, protein and proteoform nodes

### 2. Characterization of single network topology:

Metrics will be shown in tables with a color coding scheme:

- Size: Number of nodes, links, accessioned entities, small molecules
- Node degree: Cases when proteoforms have higher or lower degree with examples
- Count how often proteoforms are identical, intermediate and distinct
- Observability (percolation analysis): Quantify predicting power with threshold proteoform interactome
- Connectedness
- Articulation points
- Bridges
- Path length
- Average local clustering coefficient
- Cases when proteoform nodes have higher or lower degree with examples. (Fig)
- Examples when proteoforms have distinct set of neighbors, therefore influencing the 'guilty by association' process. (Fig)
- Select and plot examples of netowrks that change their topology when using proteoforms or small molecules

### 3. Characterization of pairwise network relationshis:

- Network separation index
- Jaccard index
- Overlap index
- Statistics on entities part of the overlap
- Examples of pairs that separate when using proteoforms
- Examples of pairs that share mostly modified proteins
