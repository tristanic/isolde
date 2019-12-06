# @Author: Tristan Croll <tic20>
# @Date:   26-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 29-Mar-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll

from .molobject import (
    get_proper_dihedral_mgr,
    get_chiral_mgr,
    get_ramachandran_mgr,
    get_rotamer_mgr
)

def get_chiral_restraint_mgr(model):
    '''
    Get the :class:`ChiralRestraintMgr` for the given model, creating it if it
    doesn't yet exist.

    Args:
        * model:
            - a :class:`chimerax.AtomicStructure`
    '''
    from .molobject import ChiralRestraintMgr
    for m in model.child_models():
        if isinstance(m, ChiralRestraintMgr):
            return m
    return ChiralRestraintMgr(model)

def get_proper_dihedral_restraint_mgr(model):
    '''
    Get the :class:`ProperDihedralRestraintMgr` for the given model, creating
    it if it doesn't yet exist.

    Args:
        * model:
            - a :class:`chimerax.AtomicStructure`
    '''
    from .molobject import ProperDihedralRestraintMgr
    for m in model.child_models():
        if isinstance(m, ProperDihedralRestraintMgr):
            return m
    return ProperDihedralRestraintMgr(model)

def get_rotamer_restraint_mgr(model):
    '''
    Get the :class:`RotamerRestraintMgr` for the given model, creating it if
    it doesn't yet exist.

    Args:
        * model:
            - a :class:`chimerax.AtomicStructure`
    '''
    from .molobject import RotamerRestraintMgr
    # RotamerRestraintMgr is subordinate to ProperDihedralRestraintMgr,
    # so we need to go deeper than just child_models()
    for m in model.all_models():
        if isinstance(m, RotamerRestraintMgr):
            return m
    return RotamerRestraintMgr(model)

def get_position_restraint_mgr(model):
    '''
    Get the :class:`PositionRestraintMgr` for the given model, creating it if
    it doesn't yet exist.

    Args:
        * model:
            - a :class:`chimerax.AtomicStructure`
    '''
    from .molobject import PositionRestraintMgr
    for m in model.child_models():
        if isinstance(m, PositionRestraintMgr):
            return m
    return PositionRestraintMgr(model)

def get_distance_restraint_mgr(model):
    '''
    Get the :class:`DistanceRestraintMgr` for the given model, creating it if
    it doesn't yet exist.

    Args:
        * model:
            - a :class:`chimerax.AtomicStructure`
    '''
    from .molobject import DistanceRestraintMgr
    for m in model.child_models():
        if isinstance(m, DistanceRestraintMgr):
            return m
    return DistanceRestraintMgr(model)

def get_adaptive_distance_restraint_mgr(model, name='Adaptive Distance Restraints'):
    '''
    Get a :class:`AdaptiveDistanceRestraintMgr` for the given model, creating
    it if it doesn't yet exist.

    Args:
        * model:
            - a :class:`chimerax.AtomicStructure`
    '''
    from .molobject import AdaptiveDistanceRestraintMgr
    for m in model.child_models():
        if isinstance(m, AdaptiveDistanceRestraintMgr) and m.name == name:
            return m
    return AdaptiveDistanceRestraintMgr(model, name=name)

def get_all_adaptive_distance_restraint_mgrs(model):
    '''
    Returns a list of all current adaptive distance restraint managers.

    Args:
        * model:
            - a :class:`chimerax.AtomicStructure`
    '''
    ret = []
    from .molobject import AdaptiveDistanceRestraintMgr
    for m in model.child_models():
        if isinstance (m, AdaptiveDistanceRestraintMgr):
            ret.append(m)
    return ret

def get_tuggable_atoms_mgr(model, allow_hydrogens=None):
    '''
    Get the :class:`TuggableAtomsMgr` for the given model, creating it if it
    doesn't yet exist.

    Args:
        * model:
            - a :class:`chimerax.AtomicStructure`
    '''
    from .molobject import TuggableAtomsMgr
    for m in model.child_models():
        if isinstance(m, TuggableAtomsMgr):
            if allow_hydrogens is not None:
                m.allow_hydrogens = allow_hydrogens
            return m
    if allow_hydrogens is None:
        allow_hydrogens = 'polar'
    return TuggableAtomsMgr(model, allow_hydrogens=allow_hydrogens)

def get_mdff_mgr(model, volume, create=False):
    '''
    Get the :class:`MDFFMgr` for the given model and volume, optionally
    creating it if it doesn't yet exist.

    Args:
        * model:
            - a :class:`chimerax.AtomicStructure` instance
        * volume:
            - a :class:`chimerax.map.Volume` instance. The :class:`MDFFMgr`
              will be added as a submodel to the :class:`Volume`
        * create (default=False):
            - if True, creates and returns a :class:`MDFFMgr` if one doesn't
              already exist
    '''
    from chimerax.map import Volume
    from chimerax.atomic import AtomicStructure
    if (not isinstance(model, AtomicStructure) or
        not isinstance(volume, Volume)):
        return None
    from .molobject import MDFFMgr
    for m in volume.all_models():
        if isinstance(m, MDFFMgr):
            if m.model == model:
                return m
    if create:
        return MDFFMgr(model, volume)
    else:
        return None

def get_rota_annotator(model, create=True):
    '''
    Get the :class:`RotamerAnnotator` for the given model, optionally creating
    it if it doesn't yet exist.

    Args:
        * model:
            - a :class:`chimerax.AtomicStructure` instance
        * create (default=True):
            - if True and no :class:`RotamerAnnotator` already exists for the
              model, one will be created and returned.
    '''
    from .validation.rota_annotation import RotamerAnnotator
    for m in model.child_models():
        if isinstance(m, RotamerAnnotator):
            return m
    if create:
        return RotamerAnnotator(model)

def get_RamaAnnotator(model, create=True):
    '''
    Get the :class:`RamaAnnotator` for the given model, optionally creating it
    if it doesn't yet exist.

    Args:
        * model:
            - a :class:`chimerax.AtomicStructure` instance
        * create (default=True):
            - if True and no :class:`RamaAnnotator` already exists for the
            - model, one will be created and returned.
    '''
    from .validation.rama_annotation import RamaAnnotator
    for m in model.child_models():
        if isinstance(m, RamaAnnotator):
            return m
    if create:
        return RamaAnnotator(model)
