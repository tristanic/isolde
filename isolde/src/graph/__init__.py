# @Author: Tristan Croll <tic20>
# @Date:   17-Jun-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 26-Jun-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

from ..mcsplit import Graph

def make_graph_from_residue(residue):
    '''
    Create a graph representation of a residue's topology, with nodes labelled
    by atomic number.
    '''
    import numpy
    from chimerax.atomic import Atoms
    ratoms = residue.atoms
    labels = ratoms.elements.numbers
    edges = numpy.array([ratoms.indices(Atoms(b.atoms)) for b in ratoms.intra_bonds])
    return Graph(labels, edges)

def make_graph_from_residue_template(tmpl):
    '''
    Create a graph representation of a TmplResidue's topology, with nodes labelled
    by atomic number.
    '''
    import numpy
    tatoms = tmpl.atoms
    labels = numpy.array([a.element.number for a in tatoms])
    bonds = set()
    for a in tmpl.atoms:
        for n in a.neighbors:
            bonds.add(frozenset([tatoms.index(aa) for aa in (a, n)]))
    edges = numpy.array([list(b) for b in bonds])
    return Graph(labels, edges)

def save_residue_graphs_to_text(residues, filename):
    import numpy
    from chimerax.atomic import Atoms
    with open(filename, 'wt') as outfile:
        for r in residues:
            ratoms = r.atoms
            labels = ratoms.elements.numbers
            edges = numpy.array([ratoms.indices(Atoms(b.atoms)) for b in ratoms.intra_bonds])
            outfile.write('graph_{}\n'.format(r.name))
            outfile.write('labels\n')
            outfile.write(', '.join([str(l) for l in labels])+'\n')
            outfile.write('edges\n')
            for e in edges:
                outfile.write('{},{}\n'.format(e[0],e[1]))
