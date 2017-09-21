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
from .dihedrals import Dihedrals

class IsoldeModel(Model):
    '''
    Base class. Should not be instantiated directly.
    '''
    def __init__(self, session, isolde, model):
        self.session = session
        self.isolde = isolde
        self.add([model])
        d = self.annotations = Drawing('ISOLDE Annotations')
        self.add_drawing(d)
        self._all_annotated_dihedrals = Dihedrals(drawing=d, session=session)

    @property
    def master_model(self):
        return self._master_model

class IsoldeCrystalModel(IsoldeModel):
    '''
    Prepares a crystal structure for ISOLDE
    '''
    def __init__(self, session, isolde, crystal_structure):
        if not is_instance(crystal_structure, CrystalStructure):
            raise TypeError('crystal_structure should be a Clipper CrystalStructure!')
        super().__init__(session, isolde, crystal_structure)
        self._master_model = crystal_structure.master_model

class IsoldeEMModel(IsoldeModel):
    '''
    Prepares an EM structure and any associated real-space maps for ISOLDE
    '''
    def __init__(self, session, isolde, atomic_structure, maps):
        super().__init__(session, isolde, atomic_structure)
        self._master_model = atomic_structure
        #TODO: create a model using the Clipper infrastructure using real-space
        #      maps

class IsoldeFreeModel(IsoldeModel):
    '''
    Prepares a molecule for simulation in ISOLDE in the absence of any maps.
    '''
    def __init__(self, session, isolde, atomic_structure):
        super().__init__(session, isolde, atomic_structure)
        self._master_model = atomic_structure
