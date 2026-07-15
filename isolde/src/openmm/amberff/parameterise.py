
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
        try:
            amber_to_ffxml(frcmod_file, ante_out, output_name=residue.name+'.xml', atom_names_from=residue)
        except UserError:
            session.logger.warning("Failed to assign atom names in the template based on the residue. This shouldn't "
                                   "affect most use cases, but will make hand-editing of the template more difficult.")
            amber_to_ffxml(frcmod_file, ante_out, output_name=residue.name+'.xml', atom_names_from=None)
        import shutil
        shutil.copyfile(ante_out, f'{residue.name}.mol2')

def parameterise_cmd(session, residues, override=False, net_charge=None,
                     always_raise_errors=True, shell_radius=1, base_templates=None):
    from chimerax.core.errors import UserError
    from chimerax.atomic import Residues

    # Metal-involved residues (a residue that contains, or is bonded to, a metal
    # atom) route to the bonded metal-site pipeline first -- GAFF2 cannot type a
    # metal, so these must not fall through to the covalent/free organic paths.
    def _metal_involved(r):
        if any(a.element.is_metal for a in r.atoms):
            return True
        return any(nb.element.is_metal for a in r.atoms for nb in a.neighbors)
    metal_seeds = [r for r in residues if _metal_involved(r)]
    handled_metal = set()
    if metal_seeds:
        from .covalent import detect_metal_site, parameterise_metal_site
        for r in metal_seeds:
            if r in handled_metal:
                continue
            try:
                site = detect_metal_site(r)
            except UserError as e:
                if always_raise_errors:
                    raise
                session.logger.warning(str(e))
                continue
            handled_metal.update(site.residues)
            try:
                parameterise_metal_site(session, site, shell_radius=shell_radius,
                                        net_charge=net_charge,
                                        base_templates=base_templates)
            except Exception as e:
                if always_raise_errors:
                    raise UserError(str(e))
                session.logger.warning('Metal-site parameterisation of %r failed: %s'
                                       % (site, e))
        residues = Residues([r for r in residues if r not in handled_metal])
        if not len(residues):
            return

    # Covalent residues (bonded to another residue) route to the covalent-unit
    # pipeline; free ligands take the classic single-residue path below.
    covalent = [r for r in residues if len(r.neighbors) != 0]
    if covalent:
        from .covalent import detect_covalent_unit, parameterise_covalent_unit
        handled = set()
        for r in covalent:
            if r in handled:
                continue
            try:
                unit = detect_covalent_unit(r)
            except UserError as e:
                if always_raise_errors:
                    raise
                session.logger.warning(str(e))
                continue
            handled.update(unit.residues)
            try:
                parameterise_covalent_unit(session, unit, shell_radius=shell_radius,
                                           net_charge=net_charge,
                                           base_templates=base_templates)
            except Exception as e:
                if always_raise_errors:
                    raise UserError(str(e))
                session.logger.warning('Covalent parameterisation of %r failed: %s'
                                       % (unit, e))

    free = Residues([r for r in residues if len(r.neighbors) == 0])
    if not len(free):
        return
    unique_residue_types = [free[free.names == name][0] for name in free.unique_names]
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
            # Free ligands now go through the RDKit pipeline (correct bond orders,
            # order-based atom-name round-trip, self-contained collision-proof
            # template) -- the same machinery as the covalent path.
            from .covalent import parameterise_free_ligand
            parameterise_free_ligand(session, residue, net_charge=net_charge,
                                     base_templates=base_templates)
        except Exception as e:
            if always_raise_errors:
                raise UserError(str(e))
            else:
                session.logger.warning(f'Parameterisation of {residue.name} failed with the following message:')
                session.logger.warning(str(e))
                continue

def register_isolde_param(logger):
    from chimerax.atomic import ResiduesArg
    from chimerax.core.commands import (CmdDesc, BoolArg, IntArg, StringArg, ListOf,
                                        register)
    desc = CmdDesc(
        required=[('residues', ResiduesArg)],
        keyword=[
            ('override', BoolArg),
            ('net_charge', IntArg),
            ('always_raise_errors', BoolArg),
            ('shell_radius', IntArg),
            ('base_templates', ListOf(StringArg)),
        ],
        synopsis=('Parameterise ligand(s) for use in ISOLDE. Supports most organic species. '
            'Covalent ligands, residue-residue crosslinks and main-chain modifications are '
            'built as a capped super-residue (select any residue of the unit). Metal sites '
            '(hemes, Zn/Mg/Mn/Ca/... coordination) are built as a bonded metal template with '
            'soft empirical coordination terms (select the metalloligand or metal). Hydrogens '
            'must be present and correct. baseTemplates takes a list of CCD ids, registered '
            'template names and/or SMILES strings whose chemistry is matched onto the '
            'ligand to fix bond orders/charges when perception is unreliable.')
    )
    register('isolde parameterise', desc, parameterise_cmd, logger=logger)
