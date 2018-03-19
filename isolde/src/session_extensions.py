
def get_proper_dihedral_mgr(session):
    if hasattr(session, 'proper_dihedral_mgr') and not session.proper_dihedral_mgr.deleted:
        return session.proper_dihedral_mgr
    from .molobject import Proper_Dihedral_Mgr
    return Proper_Dihedral_Mgr(session)

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

def get_proper_dihedral_restraint_mgr(model):
    from .molobject import Proper_Dihedral_Restraint_Mgr
    for m in model.child_models():
        if isinstance(m, Proper_Dihedral_Restraint_Mgr):
            return m
    return Proper_Dihedral_Restraint_Mgr(model)

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

def get_mdff_mgr(model, volume):
    from .molobject import MDFF_Mgr
    for m in model.all_models():
        if isinstance(m, MDFF_Mgr):
            if m.volume == volume:
                return m
    return MDFF_Mgr(model, volume)

def get_rota_annotator(model):
    from .validation_new.rota_annotation import Rotamer_Annotator
    for m in model.child_models():
        if isinstance(m, Rotamer_Annotator):
            return m
    return Rotamer_Annotator(model)

def get_rama_annotator(model):
    from .validation_new.rama_annotation import Rama_Annotator
    for m in model.child_models():
        if isinstance(m, Rama_Annotator):
            return m
    return Rama_Annotator(model)
