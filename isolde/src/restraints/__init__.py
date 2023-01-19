# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 01-Apr-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



from ..molobject import RestraintChangeTracker

# Manager classes
from ..molobject import (
    MDFFMgr,
    PositionRestraintMgr,
    TuggableAtomsMgr,
    DistanceRestraintMgr,
    AdaptiveDistanceRestraintMgr,
    ProperDihedralRestraintMgr,
    AdaptiveDihedralRestraintMgr,
    RotamerRestraintMgr,
    ChiralRestraintMgr
    )
#Singular classes
from ..molobject import (
    MDFFAtom,
    PositionRestraint,
    TuggableAtom,
    DistanceRestraint,
    AdaptiveDistanceRestraint,
    ProperDihedralRestraint,
    AdaptiveDihedralRestraint,
    RotamerRestraint,
    ChiralRestraint
    )
# Array classes
from ..molarray import (
    MDFFAtoms,
    PositionRestraints,
    TuggableAtoms,
    DistanceRestraints,
    AdaptiveDistanceRestraints,
    ProperDihedralRestraints,
    AdaptiveDihedralRestraints,
    RotamerRestraints,
    ChiralRestraints
    )
# Utility methods
from .restraint_utils import (
    restrain_torsions_to_template,
    restrain_ca_distances_to_template,
    restrain_atom_distances_to_template,
    restrain_small_ligands,
    restrain_secondary_structure
)
