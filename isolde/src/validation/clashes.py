
import numpy

SEVERE_CUTOFF=0.6
STRICT_CUTOFF=0.4

def unique_clashes(session, atoms, severe_only = False):
    from chimerax.clashes import find_clashes
    if severe_only:
        cutoff = SEVERE_CUTOFF
    else:
        cutoff = STRICT_CUTOFF
    clash_dict = find_clashes(session, atoms, inter_model=False, inter_submodel=False, clash_threshold=cutoff)
    unique_clashes = set()
    for a, cdict in clash_dict.items():
        for a2, overlap in cdict.items():
            unique_clashes.add(Clash(a, a2, overlap))
    unique_clashes = sorted(unique_clashes, key=lambda c:c.overlap, reverse=True)
    return Clashes(list(unique_clashes))


class Clash:
    def __init__(self, atom1, atom2, overlap):
        self._atoms = tuple(sorted([atom1, atom2]))
        self.overlap = overlap
    
    @property
    def atoms(self):
        from chimerax.atomic import Atoms
        return Atoms(self._atoms)

    def __hash__(self):
        return hash(self._atoms)
    
    def __eq__(self, other):
        return self._atoms == other._atoms
    
    @property
    def center(self):
        a1, a2 = self.atoms
        return (a1.coord + a2.coord)/2

class Clashes:
    def __init__(self, clashes):
        self.clashes = numpy.array(clashes)

    def __len__(self):
        return len(self.clashes)
    
    def __getitem__(self, val):
        if isinstance(val, int):
            return self.clashes[val]
        return Clashes(self.clashes[val])
    
    def __iter__(self):
        return iter(self.clashes)
    
    @property
    def centers(self):
        return numpy.array(c.center for c in self.clashes)