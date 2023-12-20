
supported_elements = ('C','N','O','S','P','H','F','Cl','Br','I')

def residue_qualifies(residue):
    unsupported_elements = set(residue.atoms.element_names).difference(set(supported_elements))
    return len(unsupported_elements==0)

def parameterise_ligand(session, residue, net_charge=None, charge_method='am1-bcc'):
    from chimerax.core.errors import UserError
    if len(residue.neighbors) != 0:
        raise UserError(f"Residue is covalently bound to {len(residue.neighbors)} other residues. "
            "ISOLDE currenly only supports parameterising new non-covalent ligands.")
    import numpy
    element_check = numpy.isin(residue.atoms.element_names, supported_elements)
    if not numpy.all(element_check):
        unsupported = residue.atoms[numpy.logical_not(element_check)]
        raise UserError('Automatic ligand parameterisation currently only supports the following elements:\n'
            f'{", ".join(supported_elements)}\n'
            f'Residue type {residue.name} contains the unsupported elements {", ".join(numpy.unique(unsupported.element_names))}.'
            'If you wish to work with this residue you will need to parameterise it using an external package.')
    import os
    base_path = os.path.dirname(os.path.abspath(__file__))
    gaff2_parms = os.path.join(base_path, 'gaff2.dat')
    from chimerax.add_charge.charge import nonstd_charge, estimate_net_charge
    from chimerax.atomic import Residues
    if net_charge is None:
        # Workaround - trigger recalculation of idatm_types if necessary to make sure rings are current
        atom_types = residue.atoms.idatm_types
        net_charge = estimate_net_charge(residue.atoms)
    from chimerax.amber_info import amber_bin
    import tempfile
    with tempfile.TemporaryDirectory() as tempdir:
        nonstd_charge(session, Residues([residue]), net_charge, charge_method, temp_dir=tempdir)
        ante_out = os.path.join(tempdir, "ante.out.mol2")
        parmchk2 = os.path.join(amber_bin, 'parmchk2')
        frcmod_file = os.path.join(tempdir, "residue.frcmod")
        import subprocess
        subprocess.check_output([parmchk2, '-s','2','-i',ante_out,'-f','mol2','-p',gaff2_parms,'-o',frcmod_file])
        from .amber_convert import amber_to_ffxml
        amber_to_ffxml(frcmod_file, ante_out, output_name=residue.name+'.xml')

def parameterise_cmd(session, residues, override=False, net_charge=None, always_raise_errors=True):
    from chimerax.core.errors import UserError
    unique_residue_types = [residues[residues.names==name][0] for name in residues.unique_names]
    for residue in unique_residue_types:
        if hasattr(session, 'isolde'):
            ff_name = session.isolde.sim_params.forcefield
            forcefield = session.isolde.forcefield_mgr[ff_name]
            ligand_db = session.isolde.forcefield_mgr.ligand_db(ff_name)
            from chimerax.isolde.openmm.openmm_interface import find_residue_templates
            from chimerax.atomic import Residues
            templates = find_residue_templates(Residues([residue]), forcefield, ligand_db=ligand_db, logger=session.logger)
            if len(templates):
                if not override:
                    raise UserError(f'Residue name {residue.name} already corresponds to template {templates[0]} in '
                        f'the {ff_name} forcefield. If you wish to replace that template, re-run this '
                        'command with override=True')
        try:
            parameterise_ligand(session, residue, net_charge=net_charge)
        except Exception as e:
            if always_raise_errors:
                raise UserError(str(e))
            else:
                session.logger.warning(f'Parameterisation of {residue.name} failed with the following message:')
                session.logger.warning(str(e))
                continue
        session.logger.info(f'OpenMM ffXML file {residue.name} written to the current working directory.')
        if hasattr(session, 'isolde'):
            if len(templates):
                forcefield._templates.pop(templates[0])
            forcefield.loadFile(f'{residue.name}.xml', resname_prefix='USER_')
            session.logger.info(f'New template added to forcefield as USER_{residue.name}. This ligand should '
                'now work in all remaining simulations for this session. To use in '
                'future sessions, load the ffXML file with ISOLDE\'s Load Residue MD Definition(s) button.')

def register_isolde_param(logger):
    from chimerax.atomic import ResiduesArg
    from chimerax.core.commands import CmdDesc, BoolArg, IntArg, register 
    desc = CmdDesc(
        required=[('residues', ResiduesArg)],
        keyword=[
            ('override', BoolArg),
            ('net_charge', IntArg),
            ('always_raise_errors', BoolArg)
        ],
        synopsis=('Parameterise ligand(s) for use in ISOLDE. Supports most organic species; '
            'covalent bonds between residues are not supported. Hydrogens must be present and correct.')
    )
    register('isolde parameterise', desc, parameterise_cmd, logger=logger)
