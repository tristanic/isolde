

def nmr_restraints_from_talos_file(session, model, chain_id, filename, spring_constant = 100):
    KAPPA_MAX=30
    from chimerax.isolde import session_extensions as sx
    from math import radians
    adrm = sx.get_adaptive_dihedral_restraint_mgr(model)
    from chimerax.atomic import Residues
    import numpy
    residues = model.residues[model.residues.chain_ids==chain_id]
    with open(filename, 'rt') as f:
        for line in f:
            comment = False
            for fw in ('REMARK', 'DATA', 'VARS', 'FORMAT'):
                if line.startswith(fw):
                    comment = True
            if comment:
                continue
            if not len(line.strip()):
                continue
            resnum = int(line[:5])
            phi = float(line[7:15])
            psi = float(line[16:24])
            dphi = float(line[25:33])
            dpsi = float(line[34:42])
            note = line[64:].strip().lower()
            if note not in ('strong','generous'):
                continue
            try:
                res = residues[residues.numbers==resnum][0]
            except IndexError:
                session.logger.warning(f'Restraint specified on residue /{chain_id}:{resnum}, but that residue does not exist!')
                continue
            phi_r = adrm.add_restraint_by_residue_and_name(res, 'phi')
            if phi_r is None:
                session.logger.warning(f'Failed to add phi restraint for residue {resnum}!')
            if phi_r is not None:
                phi_r.target = radians(phi)
                phi_r.kappa = min(radians(dphi)**(-2),KAPPA_MAX)
                phi_r.alpha = 1
                phi_r.enabled = True
            psi_r = adrm.add_restraint_by_residue_and_name(res, 'psi')
            if psi_r is not None:
                psi_r.target = radians(psi)
                psi_r.kappa = min(radians(dpsi)**(-2),KAPPA_MAX)
                psi_r.alpha = 1
                psi_r.enabled = True

def nmr_restraints_cmd(session, chain, filename, spring_constant=None):
    if len(chain) != 1:
        from chimerax.core.errors import UserError
        raise UserError('Must specify a single chain!')
    chain_id = chain[0].chain_id
    model = chain[0].structure
    if spring_constant is None:
        spring_constant = 100
    nmr_restraints_from_talos_file(session, model, chain_id, filename, spring_constant)

def register_restrain_talos_nmr(logger):
    from chimerax.core.commands import CmdDesc, register, FileNameArg, FloatArg            
    from chimerax.atomic import UniqueChainsArg

    desc = CmdDesc(
        required = [
            ('chain', UniqueChainsArg),
            ('filename', FileNameArg)
        ],
        keyword = [
            ('spring_constant', FloatArg)
        ],
        synopsis="(EXPERIMENTAL) Add torsion restraints to a chain from a TALOS-N NMR restraints file"
    )

    register ('isolde restrain talosNMR', desc, nmr_restraints_cmd, logger=logger)
