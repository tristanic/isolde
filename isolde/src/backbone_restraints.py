# Copyright 2017 Tristan Croll
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


'''
Classes and functions for handling of secondary structure restraints,
register shifts, and other operations on the protein backbone.
'''
import numpy
from chimerax.core.atomic import Residue, Atoms, AtomicStructure
from .restraints_base import Distance_Restraint, Distance_Restraints
from .restraints_base import Position_Restraint, Position_Restraints

MAX_RESTRAINT_FORCE = 250.0 # kJ/mol/A

class CA_to_CA_plus_Two(Distance_Restraints):
    '''
    Restraining the CA - (CA+2) distance is very useful when forcing
    specific secondary structures - in particular, beta strands, where
    simply forcing the backbone dihedrals is almost never enough.
    '''
    HELIX_DISTANCE = 5.43
    STRAND_DISTANCE = 6.81
    def __init__(self, session, model_or_list, name = 'ca_to_ca_plus_two'):
        '''
        Finds all connected CA - (CA+2) atom pairs in an atomic model, and
        creates a Distance_Restraints object covering all of them.
        '''
        self.session = session
        if isinstance(model_or_list, AtomicStructure):
            model = model_or_list
            pbg = self._pseudobondgroup = session.pb_manager.get_group(model.name + 'restraint pseudobonds')
            model = model_or_list
            fragments = model.polymers(
                missing_structure_treatment = model.PMS_NEVER_CONNECTS)
            restraints = []
            for f in fragments:
                if f[1] != Residue.PT_AMINO:
                    continue
                f = f[0]
                ca_atoms = f.atoms.filter(f.atoms.names == 'CA')
                for c1, c2 in zip(ca_atoms[0:-2], ca_atoms[2:]):
                    restraints.append(Distance_Restraint(Atoms((c1, c2)), name, 0, 0, pseudobond_group = pbg))
        else:
            restraints = model_or_list
            pbg = None
        super().__init__(session, restraints, name = name)

    def __getitem__(self, i):
        if isinstance(i, Residue):
            a = i.atoms
            return self[a[a.names == 'CA']]
        return super().__getitem__(i)



    pass

class O_to_N_plus_Four(Distance_Restraints):
    '''
    For restraining alpha helices, the standard approach is to apply
    distance restraints to the (n, n+4) H-bonds as well as the standard
    phi/psi backbone restraints. Better to avoid applying restraints to
    hydrogens, so here we'll connect the carbonyl oxygens to the (n+4)
    nitrogens.
    '''
    HELIX_DISTANCE = 3.05
    STRAND_DISTANCE = 11.4
    def __init__(self, session, model_or_list, name = 'o_to_n_plus_four'):
        '''
        Finds all connected O - (N+4) atom pairs in an atomic model, and
        creates a Distance_Restraints object covering all of them.
        '''
        self.session=session
        if isinstance(model_or_list, AtomicStructure):
            model = model_or_list
            pbg = self._pseudobondgroup = session.pb_manager.get_group(model.name + 'restraint pseudobonds')
            fragments = model.polymers(
                missing_structure_treatment = model.PMS_NEVER_CONNECTS)
            restraints = []
            for f in fragments:
                if f[1] != Residue.PT_AMINO:
                    continue
                f = f[0]
                o_atoms = f.atoms.filter(f.atoms.names == 'O')
                n_atoms = f.atoms.filter(f.atoms.names == 'N')
                for o, n in zip(o_atoms[0:-4], n_atoms[4:]):
                    restraints.append(Distance_Restraint(Atoms((o, n)), name, 0, 0, pseudobond_group = pbg))
        else:
            restraints = model_or_list
        super().__init__(session, restraints, name = name)

    def __getitem__(self, i):
        if isinstance(i, Residue):
            a = i.atoms
            return self[a[numpy.in1d(a.names, ['N','O'])]]
        return super().__getitem__(i)
