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

        self._proper_dihedral_mgr = None
        self._rama_mgr = None
        self._rota_mgr = None

        self._rama_annotator = None
        self._rota_annotator = None

        self._proper_dihedral_restraint_mgr = None
        self._position_restraint_mgr = None
        self._tuggable_atoms_mgr = None
        self._distance_restraint_mgr = None

        session.models.add([self])


    @property
    def master_model(self):
        return self._master_model

    @property
    def proper_dihedral_mgr(self):
        if self._proper_dihedral_mgr is None:
            from . import session_extensions
            self._proper_dihedral_mgr = session_extensions.get_proper_dihedral_mgr(self.session)
        return self._proper_dihedral_mgr

    @property
    def rama_mgr(self):
        if self._rama_mgr is None:
            from . import session_extensions
            self._rama_mgr = session_extensions.get_ramachandran_mgr(self.session)
        return self._rama_mgr

    @property
    def rota_mgr(self):
        if self._rota_mgr is None:
            from . import session_extensions
            self._rota_mgr = session_extensions.get_rotamer_mgr(self.session)

    @property
    def rama_annotator(self):
        if self._rama_annotator is None:
            from . import session_extensions
            self._rama_annotator = session_extensions.get_rama_annotator(self.master_model)
        return self._rama_annotator

    @property
    def rota_annotator(self):
        if self._rota_annotator is None:
            from . import session_extensions
            self._rota_annotator = session_extensions.get_rota_annotator(self.master_model)
        return self._rota_annotator

    @property
    def proper_dihedral_restraint_mgr(self):
        if self._proper_dihedral_restraint_mgr is None:
            from . import session_extensions
            self._proper_dihedral_restraint_mgr = session_extensions.get_proper_dihedral_restraint_mgr(self.master_model)
        return self._proper_dihedral_restraint_mgr

    @property
    def position_restraint_mgr(self):
        if self._position_restraint_mgr is None:
            from . import session_extensions
            self._position_restraint_mgr = session_extensions.get_position_restraint_mgr(self.master_model)
        return self._position_restraint_mgr

    @property
    def tuggable_atoms_mgr(self):
        if self._tuggable_atoms_mgr is None:
            from . import session_extensions
            self._tuggable_atoms_mgr = session_extensions.get_tuggable_atoms_mgr(self.master_model)
        return self._tuggable_atoms_mgr

    @property
    def distance_restraint_mgr(self):
        if self._distance_restraint_mgr is None:
            from . import session_extensions
            self._distance_restraint_mgr = session_extensions.get_distance_restraint_mgr(self.master_model)
        return self._distance_restraint_mgr


def _isolde_crystal_model_from_atomic_structure_and_mtz(model, mtzfile,
                                        map_oversampling=3):
    session = model.session
    cs = CrystalStructure(session, model, mtzfile=mtzfile, map_oversampling = map_oversampling)
    return IsoldeCrystalModel(cs)

class IsoldeCrystalModel(IsoldeModel):
    '''
    Prepares a crystal structure for ISOLDE
    '''
    def __init__(self, atomic_structure, mtz_file, map_oversampling = 3.0):
        self._master_model = atomic_structure
        from chimerax.clipper import symmetry
        self._map_handler = symmetry.XtalSymmetryHandler(atomic_structure, mtz_file,
            map_oversampling=map_oversampling)
        super().__init__(atomic_structure)

    @property
    def map_handler(self):
        return self._map_handler

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
