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
    def __init__(self, model):
        if self.__class__ == IsoldeModel:
            raise RuntimeError('IsoldeModel is a base class and should not be '\
                +'instantiated directly. Use one of the derived classes: '\
                +'IsoldeCrystalModel, IsoldeEMModel or IsoldeFreeModel.')
        session = model.session
        super().__init__('ISOLDE ' + model.name, session)
        self.add([model])
        # _master_model is defined by the derived class
        mm = self._master_model

        # Initialise validation managers
        from .validation_new import Rama_Annotator, Rotamer_Annotator
        rama_a = self._rama_annotator = Rama_Annotator(session, mm)
        rota_a = self._rota_annotator = Rotamer_Annotator(session, mm)

        # Initialise restraint managers
        from .molobject import Position_Restraint_Mgr, Distance_Restraint_Mgr, \
                            Proper_Dihedral_Restraint_Mgr
        pr_mgr = self._position_restraint_mgr = Position_Restraint_Mgr(mm)
        distr_mgr = self._distance_restraint_mgr = Distance_Restraint_Mgr(mm)
        pdr_mgr = self._proper_dihedral_restraint_mgr = Proper_Dihedral_Restraint_Mgr(mm)


        session.models.add([self])


    @property
    def master_model(self):
        return self._master_model

    @property
    def rama_annotator(self):
        return self._rama_annotator

    @property
    def rota_annotator(self):
        return self_rota_annotator

    @property
    def position_restraint_mgr(self):
        return self._position_restraint_mgr

    @property
    def distance_restraint_mgr(self):
        return self._distance_restraint_mgr

    @property
    def proper_dihedral_restraint_mgr(self):
        return self._proper_dihedral_restraint_mgr



def _isolde_crystal_model_from_atomic_structure_and_mtz(model, mtzfile,
                                        map_oversampling=3):
    session = model.session
    cs = CrystalStructure(session, model, mtzfile=mtzfile, map_oversampling = map_oversampling)
    return IsoldeCrystalModel(cs)

class IsoldeCrystalModel(IsoldeModel):
    '''
    Prepares a crystal structure for ISOLDE
    '''
    def __init__(self, crystal_structure):
        if not isinstance(crystal_structure, CrystalStructure):
            raise TypeError('crystal_structure should be a Clipper CrystalStructure!')
        self._master_model = crystal_structure.master_model
        super().__init__(crystal_structure)

class IsoldeEMModel(IsoldeModel):
    '''
    Prepares an EM structure and any associated real-space maps for ISOLDE
    '''
    def __init__(self, atomic_structure, maps):
        self._master_model = atomic_structure
        super().__init__(atomic_structure)
        #TODO: create a model using the Clipper infrastructure using real-space
        #      maps

class IsoldeFreeModel(IsoldeModel):
    '''
    Prepares a molecule for simulation in ISOLDE in the absence of any maps.
    '''
    def __init__(self, atomic_structure):
        self._master_model = atomic_structure
        super().__init__(atomic_structure)
