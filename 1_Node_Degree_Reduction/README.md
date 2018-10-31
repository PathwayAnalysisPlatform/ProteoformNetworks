# Data

As part of the arguments to prove that the proteoforms networks are better suited for functional
analysis we analyse the degree of the nodes.

Interaction networks are composed of a set of nodes and links.
For protein interaction networks, the nodes represent an protein accessions and links are interactions among two proteins.
For the proteoforms interaction networks, there is one node for each annotated proteoform of a protein. Since isoform variants and
PTMs change the structure and function of the protein the proteoforms may link to only a subset of the interacting partners of
the generic protein. This reduces the number of links each proteoform node has.

Overall, the biological networks have a share of nodes with high degree connecting to plenty of other entities, creating a 
small world effect, that reduces the length of the shortest paths in the overall network and reduces the chances to create hypothesis
based solely on the network topology. The smaller degree of proteoform nodes reduces the impact of the small world effect, improving
the accuracy of the network analysis algorithms.

* The network files are located at: PathwayMatcher\resources\networks\all\1.8.1
