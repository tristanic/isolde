# @Author: Tristan Croll <tic20>
# @Date:   26-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 26-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll



def get_proper_dihedral_mgr(session):
    '''
    Get the session-level :class:`Proper_Dihedral_Mgr` singleton, creating it if
    it doesn't yet exist.

    Args;
        * session:
            - the top-level ChimeraX session instance
    '''
    if hasattr(session, 'proper_dihedral_mgr') and not session.proper_dihedral_mgr.deleted:
        return session.proper_dihedral_mgr
    from .molobject import Proper_Dihedral_Mgr
    return Proper_Dihedral_Mgr(session)

def get_chiral_mgr(session):
    '''
    Get the session-level :class:`Chiral_Mgr` singleton, creating it if it
    doesn't yet exist.

    Args;
        * session:
            - the top-level ChimeraX session instance
    '''
    if hasattr(session, 'chiral_mgr') and not session.chiral_mgr.deleted:
        return session.chiral_mgr
    from .molobject import Chiral_Mgr
    return Chiral_Mgr(session)

def get_ramachandran_mgr(session):
    '''
    Get the session-level :class:`Rama_Mgr` singleton, creating it if it doesn't
    yet exist.

    Args;
        * session:
            - the top-level ChimeraX session instance
    '''
    if hasattr(session, 'rama_mgr') and not session.rama_mgr.deleted:
        return session.rama_mgr
    from .molobject import Rama_Mgr
    return Rama_Mgr(session)

def get_rotamer_mgr(session):
    '''
    Get the session-level :class:`Rota_Mgr` singleton, creating it if it doesn't
    yet exist.

    Args;
        * session:
            - the top-level ChimeraX session instance
    '''
    if hasattr(session, 'rota_mgr') and not session.rota_mgr.deleted:
        return session.rota_mgr
    from .molobject import Rota_Mgr
    return Rota_Mgr(session)

def get_chiral_restraint_mgr(model):
    '''
    Get the :class:`Chiral_Restraint_Mgr` for the given model, creating it if it
    doesn't yet exist.

    Args;
        * model:
            - a :class:`chimerax.AtomicStructure`
    '''
    from .molobject import Chiral_Restraint_Mgr
    for m in model.child_models():
        if isinstance(m, Chiral_Restraint_Mgr):
            return m
    return Chiral_Restraint_Mgr(model)

def get_proper_dihedral_restraint_mgr(model):
    '''
    Get the :class:`Proper_Dihedral_Restraint_Mgr` for the given model, creating
    it if it doesn't yet exist.

    Args;
        * model:
            - a :class:`chimerax.AtomicStructure`
    '''
    from .molobject import Proper_Dihedral_Restraint_Mgr
    for m in model.child_models():
        if isinstance(m, Proper_Dihedral_Restraint_Mgr):
            return m
    return Proper_Dihedral_Restraint_Mgr(model)

def get_rotamer_restraint_mgr(model):
    '''
    Get the :class:`Rotamer_Restraint_Mgr` for the given model, creating it if
    it doesn't yet exist.

    Args;
        * model:
            - a :class:`chimerax.AtomicStructure`
    '''
    from .molobject import Rotamer_Restraint_Mgr
    # Rotamer_Restraint_Mgr is subordinate to Proper_Dihedral_Restraint_Mgr,
    # so we need to go deeper than just child_models()
    for m in model.all_models():
        if isinstance(m, Rotamer_Restraint_Mgr):
            return m
    return Rotamer_Restraint_Mgr(model)

def get_position_restraint_mgr(model):
    '''
    Get the :class:`Position_Restraint_Mgr` for the given model, creating it if
    it doesn't yet exist.

    Args;
        * model:
            - a :class:`chimerax.AtomicStructure`
    '''
    from .molobject import Position_Restraint_Mgr
    for m in model.child_models():
        if isinstance(m, Position_Restraint_Mgr):
            return m
    return Position_Restraint_Mgr(model)

def get_distance_restraint_mgr(model):
    '''
    Get the :class:`Distance_Restraint_Mgr` for the given model, creating it if
    it doesn't yet exist.

    Args;
        * model:
            - a :class:`chimerax.AtomicStructure`
    '''
    from .molobject import Distance_Restraint_Mgr
    for m in model.child_models():
        if isinstance(m, Distance_Restraint_Mgr):
            return m
    return Distance_Restraint_Mgr(model)

def get_tuggable_atoms_mgr(model):
    '''
    Get the :class:`Tuggable_Atoms_Mgr` for the given model, creating it if it
    doesn't yet exist.

    Args;
        * model:
            - a :class:`chimerax.AtomicStructure`
    '''
    from .molobject import Tuggable_Atoms_Mgr
    for m in model.child_models():
        if isinstance(m, Tuggable_Atoms_Mgr):
            return m
    return Tuggable_Atoms_Mgr(model)

def get_mdff_mgr(model, volume, create=False):
    '''
    Get the :class:`MDFF_Mgr` for the given model and volume, optionally
    creating it if it doesn't yet exist.

    Args;
        * model:
            - a :class:`chimerax.AtomicStructure` instance
        * volume:
            - a :class:`chimerax.map.Volume` instance. The :class:`MDFF_Mgr`
              will be added as a submodel to the :class:`Volume`
        * create (default=False):
            - if True, creates and returns a :class:`MDFF_Mgr` if one doesn't
              already exist
    '''
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

def get_rota_annotator(model, create=True):
    '''
    Get the :class:`Rotamer_Annotator` for the given model, optionally creating
    it if it doesn't yet exist.

    Args;
        * model:
            - a :class:`chimerax.AtomicStructure` instance
        * create (default=True):
            - if True and no :class:`Rotamer_Annotator` already exists for the
              model, one will be created and returned.
    '''
    from .validation.rota_annotation import Rotamer_Annotator
    for m in model.child_models():
        if isinstance(m, Rotamer_Annotator):
            return m
    if create:
        return Rotamer_Annotator(model)

def get_rama_annotator(model, create=True):
    '''
    Get the :class:`Rama_Annotator` for the given model, optionally creating it
    if it doesn't yet exist.

    Args;
        * model:
            - a :class:`chimerax.AtomicStructure` instance
        * create (default=True):
            - if True and no :class:`Rama_Annotator` already exists for the
            - model, one will be created and returned.
    '''
    from .validation.rama_annotation import Rama_Annotator
    for m in model.child_models():
        if isinstance(m, Rama_Annotator):
            return m
    if create:
        return Rama_Annotator(model)
