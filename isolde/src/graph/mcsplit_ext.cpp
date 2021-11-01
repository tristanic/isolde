/**
 * @Author: Tristan Croll <tic20>
 * @Date:   17-Jun-2020
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 06-Jul-2020
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2016-2019 Tristan Croll
 */

// #if defined(_MSC_VER) && _MSC_VER >= 1600
// #include <string>
// const std::basic_string<char>::size_type std::basic_string<char>::npos = (std::basic_string<char>::size_type) -1;
// #endif

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "mcsp.h"

namespace py=pybind11;

PYBIND11_MODULE(mcsplit, m) {
    m.doc() = "Implementation of the McSplit algorithm for finding the maximum "
        "common induced subgraph of two graphs. https://doi.org/10.24963/ijcai.2017/99 ";

    py::class_<Graph>(m, "Graph")
        .def(py::init<unsigned int>())
        .def(py::init([](py::array_t<unsigned int> labels, py::array_t<int> edges)
        {
            if (labels.ndim() !=1)
                throw std::runtime_error("labels should be a 1-dimensional array of unsigned ints!");
            if (edges.ndim() != 2)
                throw std::runtime_error("edges should be a n x 2 array of unsigned ints!");
            std::vector<unsigned int> lab(labels.shape(0));
            std::vector<int> edges1(edges.shape(0));
            std::vector<int> edges2(edges.shape(0));
            for (size_t i=0; i<labels.shape(0); ++i)
            {
                lab[i] = labels.at(i);
            }
            for (size_t i=0; i<edges.shape(0); ++i)
            {
                edges1[i] = edges.at(i,0);
                edges2[i] = edges.at(i,1);
            }
            return make_graph(lab, edges1, edges2);
        }))
        .def_property_readonly("labels", [](const Graph& self)
        {
            return py::array(self.label.size(), self.label.data());
        })
        .def_property_readonly("adjacency_matrix", [](const Graph& self)
        {
            auto n1 = self.adjmat.size();
            if (n1==0)
                throw std::out_of_range("Adjacency matrix is empty!");
            auto n2 = self.adjmat[0].size();
            auto ret = py::array_t<unsigned int>({n1, n2});
            auto r = ret.mutable_unchecked<2>();
            for (size_t i=0; i<n1; ++i)
                for (size_t j=0; j<n2; ++j)
                    r(i,j) = self.adjmat[i][j];
            return ret;
        })
        .def_property_readonly("n", [](const Graph& self) { return self.n; })
        .def("maximum_common_subgraph", [](Graph& self, Graph& other, float timeout,
            bool connected, bool directed, bool vertex_labelled, bool edge_labelled,
            bool big_first, bool verbose, bool quiet)
        {
            std::vector<int> i0, i1;
            bool aborted;
            std::tie(i0, i1, aborted) = maximum_common_subgraph(self, other, timeout, connected,
                directed, vertex_labelled, edge_labelled, big_first, verbose, quiet);
            return py::make_tuple(py::array(i0.size(), i0.data()), py::array(i1.size(), i1.data()), aborted);
        },
        py::arg("g1"),
        py::arg("timeout")=1.0,
        py::arg("connected")=true,
        py::arg("directed")=false,
        py::arg("vertex_labelled")=true,
        py::arg("edge_labelled")=false,
        py::arg("big_first")=false,
        py::arg("verbose")=false,
        py::arg("quiet")=true
        )
        .def(py::pickle(
            [](const Graph& g) { // __getstate__
                std::vector<int> edges1, edges2;
                size_t i=0;
                // Adjacency matrix is symmetric, so we only need to iterate
                // through half of it.
                for (const auto& column: g.adjmat) {
                    for (size_t j=0; j<i; j++)
                        if (column[j] != 0)
                        {
                            edges1.push_back(i);
                            edges2.push_back(j);
                        }
                    i++;
                }
                return std::tuple<std::vector<unsigned int>, std::vector<int>, std::vector<int>>(g.label, edges1, edges2);
            },
            [](std::tuple<std::vector<unsigned int>, std::vector<int>, std::vector<int>> t) {
                std::vector<unsigned int> labels;
                std::vector<int> edges1, edges2;
                std::tie(labels, edges1, edges2) = t;
                return make_graph(labels, edges1, edges2);
                //return make_graph(t[0], t[1], t[2]);
            }
        ))
        ;

};
