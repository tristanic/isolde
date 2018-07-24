import sys, os
sys.path.insert(0, os.path.abspath('.'))
import clipper_python as clipper
import numpy

names = ['H']*5
coords = numpy.random.rand(5,3)
ua = numpy.random.rand(5,6)
occ = numpy.random.rand(5)
ui = numpy.random.rand(5)
al = clipper.Atom_list(5, names, coords, ua, occ, ui)

atom = al [2]
print(atom.element)
al.delete_atom(2)
print(atom.element)
