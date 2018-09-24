//============================================================================
// Name        : Overlap.cpp
// Author      : Luis Francisco Hernández Sánchez
// Version     : 0.0.1
// Copyright   : Your copyright notice
// Description : Read igraph network files, create a matrix with the number of
// overlapping nodes for each pair
//============================================================================

#include <iostream>
#include <igraph.h>

using namespace std;

int main(void)
{
     igraph_integer_t diameter;
     igraph_t graph;
     igraph_rng_seed(igraph_rng_default(), 42);
     igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNP, 1000, 5.0/1000,
                             IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
     igraph_diameter(&graph, &diameter, 0, 0, 0, IGRAPH_UNDIRECTED, 1);
     printf("Diameter of a random graph with average degree 5: %d\n",
             (int) diameter);
     igraph_destroy(&graph);
     return 0;
}
