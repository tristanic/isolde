#pragma once
//
// #include <stdbool.h>

#include <vector>
#include <memory>

class Graph {
public:
    Graph(unsigned int n);
    int n;
    std::vector<std::vector<unsigned int>> adjmat;
    std::vector<unsigned int> label;
};

Graph induced_subgraph(struct Graph& g, std::vector<int> vv);

std::unique_ptr<Graph> make_graph(const std::vector<unsigned int>& labels,
    const std::vector<int>& edges1,
    const std::vector<int>& edges2);
