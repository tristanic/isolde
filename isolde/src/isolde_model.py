'''
ISOLDE adds a number of annotations to any structure it considers/simulates,
to visualise restraints, cis/twisted peptide bonds, etc. In order to allow
ISOLDE to switch between models without throwing all this information out, the
necessary objects need to be bundled together with each model. Since the
ChimeraX Structure class is not designed to have child models/drawings beneath
it, this means we have to move the structure to become the child of a new model.
'''

from chimerax.core.models import Model, Drawing
from chimerax.clipper import CrystalStructure
from .dihedrals import Dihedrals, Backbone_Dihedrals
from . import position_restraints
from . import backbone_restraints
from . import rotamers


class IsoldeModel(Model):
    '''
    Base class. Should not be instantiated directly.
    '''
    def __init__(self, session, isolde, model):
        if self.__class__ == IsoldeModel:
            raise RuntimeError('IsoldeModel is a base class and should not be '\
                +'instantiated directly. Use one of the derived classes: '\
                +'IsoldeCrystalModel, IsoldeEMModel or IsoldeFreeModel.')
        super().__init__('ISOLDE ' + model.name, session)
        self.isolde = isolde
        self.add([model])
        # _master_model is defined by the derived class
        mm = self._master_model
        am = self._annotation_model = Model('ISOLDE Annotations', session)
        dm = self._dihedral_annotations_model = Model('Dihedral restraints', session)
        dd = self._dihedral_annotations_drawing = Drawing('Dihedral restraints')
        dm.add_drawing(dd)
        am.add([dm])
        self.add([am])
        ad = self._all_annotated_dihedrals = Dihedrals(drawing=dd, session=session)
        bd = self.backbone_dihedrals = Backbone_Dihedrals(session, mm)
        ad.append(bd.phi)
        ad.append(bd.psi)
        rots = self.rotamers = rotamers.all_rotamers_in_selection(session,
                                                    mm.atoms)
        for r in rots.values():
            ad.append(r.dihedrals)
        self.distance_restraints = {
            'ca_to_ca_plus_two':    backbone_restraints.CA_to_CA_plus_Two(session, mm),
            'o_to_n_plus_four':     backbone_restraints.O_to_N_plus_Four(session, mm),
        }

        self.position_restraints = position_restraints.Atom_Position_Restraints(
            session, mm, mm.atoms, triggers=isolde.triggers, create_target=True
        )


    @property
    def master_model(self):
        return self._master_model

class IsoldeCrystalModel(IsoldeModel):
    '''
    Prepares a crystal structure for ISOLDE
    '''
    def __init__(self, session, isolde, crystal_structure):
        if not isinstance(crystal_structure, CrystalStructure):
            raise TypeError('crystal_structure should be a Clipper CrystalStructure!')
        self._master_model = crystal_structure.master_model
        super().__init__(session, isolde, crystal_structure)

class IsoldeEMModel(IsoldeModel):
    '''
    Prepares an EM structure and any associated real-space maps for ISOLDE
    '''
    def __init__(self, session, isolde, atomic_structure, maps):
        self._master_model = atomic_structure
        super().__init__(session, isolde, atomic_structure)
        #TODO: create a model using the Clipper infrastructure using real-space
        #      maps

class IsoldeFreeModel(IsoldeModel):
    '''
    Prepares a molecule for simulation in ISOLDE in the absence of any maps.
    '''
    def __init__(self, session, isolde, atomic_structure):
        self._master_model = atomic_structure
        super().__init__(session, isolde, atomic_structure)
