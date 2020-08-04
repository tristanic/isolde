# @Author: Tristan Croll <tic20>
# @Date:   17-Jun-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 30-Jul-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

from ..mcsplit import Graph

def make_graph_from_residue(residue, label='element'):
    '''
    Create a graph representation of a residue's topology, with nodes labelled
    by either atomic number or an integer representation of their name.
    '''
    import numpy
    from chimerax.atomic import Atoms
    ratoms = residue.atoms
    if label=='element':
        labels = ratoms.elements.numbers
    elif label=='name':
        labels = numpy.array([simple_name_to_int(name) for name in ratoms.names])
    else:
        raise TypeError('"label" argument should be one of "element" or "name"!')
    edges = numpy.array([ratoms.indices(Atoms(b.atoms)) for b in ratoms.intra_bonds])
    return Graph(labels, edges)

def make_graph_from_residue_template(tmpl, label='element'):
    '''
    Create a graph representation of a TmplResidue's topology, with nodes labelled
    by atomic number.
    '''
    import numpy
    tatoms = tmpl.atoms
    if label=='element':
        labels = numpy.array([a.element.number for a in tatoms])
    elif label=='name':
        labels = numpy.array([simple_name_to_int(a.name) for a in tatoms])
    else:
        raise TypeError('"label" argument should be one of "element" or "name"!')
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

def save_graphs_to_text(graph_dict, filename):
    with open(filename, 'wt') as outfile:
        for graph_name, graph in graph_dict.items():
            labels, edges1, edges2 = graph.__getstate__()
            outfile.write('graph_{}\n'.format(graph_name))
            outfile.write('labels\n')
            outfile.write(', '.join([str(l) for l in labels])+'\n')
            outfile.write('edges\n')
            for e0, e1 in zip(edges1, edges2):
                outfile.write('{},{}\n'.format(e0,e1))


def simple_name_to_int(name):
    '''
    Convert an atom name to an integer for labelling graph nodes
    '''
    if len(name) > 4:
        raise RuntimeError('This method is only suitable for atom names with 4 or fewer characters!')
    multiplier = 10**((len(name)-1)*2)
    val = 0
    for char in name:
        val += ord(char)*multiplier
        multiplier /= 100
    return val
