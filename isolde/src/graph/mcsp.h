/**
 * @Author: Tristan Croll <tic20>
 * @Date:   17-Jun-2020
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 18-Jun-2020
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2016-2019 Tristan Croll
 */
#pragma once

#include "graph.h"
std::tuple<std::vector<int>, std::vector<int>, bool>
maximum_common_subgraph(Graph& g0, Graph& g1, float timeout=0,
    bool connected=true, bool directed=false, bool vertex_labelled=true,
    bool edge_labelled=true, bool big_first=false, bool verbose=false,
    bool quiet=true);
