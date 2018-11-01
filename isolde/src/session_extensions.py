# @Author: Tristan Croll <tic20>
# @Date:   26-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 26-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll



def get_proper_dihedral_mgr(session):
    if hasattr(session, 'proper_dihedral_mgr') and not session.proper_dihedral_mgr.deleted:
        return session.proper_dihedral_mgr
    from .molobject import Proper_Dihedral_Mgr
    return Proper_Dihedral_Mgr(session)

def get_chiral_mgr(session):
    if hasattr(session, 'chiral_mgr') and not session.chiral_mgr.deleted:
        return session.chiral_mgr
    from .molobject import Chiral_Mgr
    return Chiral_Mgr(session)

def get_ramachandran_mgr(session):
    if hasattr(session, 'rama_mgr') and not session.rama_mgr.deleted:
        return session.rama_mgr
    from .molobject import Rama_Mgr
    return Rama_Mgr(session)

def get_rotamer_mgr(session):
    if hasattr(session, 'rota_mgr') and not session.rota_mgr.deleted:
        return session.rota_mgr
    from .molobject import Rota_Mgr
    return Rota_Mgr(session)

def get_chiral_restraint_mgr(model):
    from .molobject import Chiral_Restraint_Mgr
    for m in model.child_models():
        if isinstance(m, Chiral_Restraint_Mgr):
            return m
    return Chiral_Restraint_Mgr(model)

def get_proper_dihedral_restraint_mgr(model):
    from .molobject import Proper_Dihedral_Restraint_Mgr
    for m in model.child_models():
        if isinstance(m, Proper_Dihedral_Restraint_Mgr):
            return m
    return Proper_Dihedral_Restraint_Mgr(model)

def get_rotamer_restraint_mgr(model):
    from .molobject import Rotamer_Restraint_Mgr
    for m in model.child_models():
        if isinstance(m, Rotamer_Restraint_Mgr):
            return m
    return Rotamer_Restraint_Mgr(model)

def get_position_restraint_mgr(model):
    from .molobject import Position_Restraint_Mgr
    for m in model.child_models():
        if isinstance(m, Position_Restraint_Mgr):
            return m
    return Position_Restraint_Mgr(model)

def get_distance_restraint_mgr(model):
    from .molobject import Distance_Restraint_Mgr
    for m in model.child_models():
        if isinstance(m, Distance_Restraint_Mgr):
            return m
    return Distance_Restraint_Mgr(model)

def get_tuggable_atoms_mgr(model):
    from .molobject import Tuggable_Atoms_Mgr
    for m in model.child_models():
        if isinstance(m, Tuggable_Atoms_Mgr):
            return m
    return Tuggable_Atoms_Mgr(model)

def get_mdff_mgr(model, volume, create=False):
    from chimerax.map import Volume
    from chimerax.atomic import AtomicStructure
    if (not isinstance(model, AtomicStructure) or
        not isinstance(volume, Volume)):
        return None
    from .molobject import MDFF_Mgr
    for m in volume.all_models():
        if isinstance(m, MDFF_Mgr):
            if m.model == model:
                return m
    if create:
        return MDFF_Mgr(model, volume)
    else:
        return None

def get_rota_annotator(model):
    from .validation.rota_annotation import Rotamer_Annotator
    for m in model.child_models():
        if isinstance(m, Rotamer_Annotator):
            return m
    return Rotamer_Annotator(model)

def get_rama_annotator(model):
    from .validation.rama_annotation import Rama_Annotator
    for m in model.child_models():
        if isinstance(m, Rama_Annotator):
            return m
    return Rama_Annotator(model)
