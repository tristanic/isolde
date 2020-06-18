# @Author: Tristan Croll <tic20>
# @Date:   17-Jun-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 17-Jun-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll
def test_make_graph(residue):
    import numpy
    from chimerax.atomic import Atoms
    ratoms = residue.atoms
    labels = ratoms.elements.numbers
    edges = numpy.array([ratoms.indices(Atoms(b.atoms)) for b in ratoms.intra_bonds])
    from ..mcsplit import Graph
    return Graph(labels, edges)
