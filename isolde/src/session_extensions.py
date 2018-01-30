
def get_proper_dihedral_manager(session):
    if hasattr(session, 'proper_dihedral_mgr') and not session.proper_dihedral_mgr.deleted:
        return session.proper_dihedral_mgr
    from .molobject import Proper_Dihedral_Mgr
    return Proper_Dihedral_Mgr(session)

def get_ramachandran_manager(session):
    if hasattr(session, 'rama_mgr') and not session.rama_mgr.deleted:
        return session.rama_mgr
    from .molobject import Rama_Mgr
    return Rama_Mgr(session)

def get_rotamer_manager(session):
    if hasattr(session, 'rota_mgr') and not session.rota_mgr.deleted:
        return session.rota_mgr
    from .molobject import Rota_Mgr
    return Rota_Mgr(session)


