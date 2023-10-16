#include "graph.h"

// #include <stdio.h>
// #include <stdlib.h>
#include <limits.h>

#include <iostream>
#include <string>
// #include <stdexcept>

constexpr int BITS_PER_UNSIGNED_INT (CHAR_BIT * sizeof(unsigned int));


Graph::Graph(unsigned int n) {
    this->n = n;
    label = std::vector<unsigned int>(n, 0u);
    adjmat = {n, std::vector<unsigned int>(n, false)};
}

Graph induced_subgraph(Graph& g, std::vector<int> vv) {
    Graph subg(vv.size());
    for (int i=0; i<subg.n; i++)
        for (int j=0; j<subg.n; j++)
            subg.adjmat[i][j] = g.adjmat[vv[i]][vv[j]];

    for (int i=0; i<subg.n; i++)
        subg.label[i] = g.label[vv[i]];
    return subg;
}

void add_edge(Graph& g, int v, int w, bool directed=false, unsigned int val=1) {
    if (v != w) {
        if (directed) {
            g.adjmat.at(v).at(w) |= val;
            g.adjmat[w][v] |= (val<<16);
        } else {
            g.adjmat.at(v).at(w) = val;
            g.adjmat[w][v] = val;
        }
    } else {
        // To indicate that a vertex has a loop, we set the most
        // significant bit of its label to 1
        g.label.at(v) |= (1u << (BITS_PER_UNSIGNED_INT-1));
    }
}

std::unique_ptr<Graph> make_graph(const std::vector<unsigned int>& labels,
    const std::vector<int>& edges1,
    const std::vector<int>& edges2)
{
    Graph* g = new Graph(labels.size());
    g->label = labels;
    for (size_t i=0; i<edges1.size(); ++i)
        add_edge(*g, edges1[i], edges2[i]);
    return std::unique_ptr<Graph>(g);
}
