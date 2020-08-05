# @Author: Tristan Croll <tic20>
# @Date:   16-Apr-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 05-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

def replace_residue(session, residue, new_residue_name):
    from chimerax.core.errors import UserError
    from chimerax.isolde.cmd import isolde_start
    isolde_start(session)
    from chimerax.atomic import Residues
    if isinstance(residue, Residues):
        if len(residue) != 1:
            raise UserError('Must have a single residue selected!')
        residue = residue[0]
    if len(residue.neighbors):
        raise UserError('Replacement by graph matching is currently only '
            'supported for ligands with no covalent bonds to other residues!')
    from ..template_utils import (
        fix_residue_from_template,
        fix_residue_to_match_md_template
    )
    from chimerax.atomic import mmcif
    try:
        cif_template = mmcif.find_template_residue(session, new_residue_name)
    except:
        err_str = ('Could not find a mmCIF template for residue name {}. '
            'For novel residues not in the Chemical Components Dictionary, '
            'you will need to provide this first.').format(new_residue_name)
        raise UserError(err_str)
    fix_residue_from_template(residue, cif_template, rename_residue=True,
        match_by='element')
    ff = session.isolde.forcefield_mgr[session.isolde.sim_params.forcefield]
    ligand_db = session.isolde.forcefield_mgr.ligand_db(session.isolde.sim_params.forcefield)
    from chimerax.isolde.openmm.openmm_interface import find_residue_templates
    from chimerax.atomic import Residues
    tdict = find_residue_templates(Residues([residue]), ff, ligand_db=ligand_db, logger=session.logger)
    md_template_name = tdict.get(0)
    if md_template_name is not None:
        fix_residue_to_match_md_template(session, residue, ff._templates[md_template_name], cif_template=cif_template)
    from chimerax.atomic import Atoms, Atom
    chiral_centers = Atoms([a for a in residue.atoms if residue.ideal_chirality(a.name) != 'N'])
    if len(chiral_centers):
        warn_str = ('Rebuilt ligand {} has chiral centres at atoms {} '
            '(highlighted). Since the current algorithm used to match topologies '
            'is not chirality aware, you should check these sites carefully to '
            'ensure they are sensible. If in doubt, it is best to delete with '
            '"del #{}/{}:{}{}" and replace with "isolde add ligand {}".').format(
                residue.name, ','.join(chiral_centers.names),
                residue.structure.id_string, residue.chain_id, residue.number,
                residue.insertion_code, residue.name)
        session.selection.clear()
        chiral_centers.selected=True
        chiral_centers.draw_modes = Atom.BALL_STYLE
        from chimerax.isolde.view import focus_on_selection
        focus_on_selection(session, chiral_centers)
        session.logger.warning(warn_str)

def add_ligand(session, *args, **kwargs):
    from chimerax.isolde.cmd import isolde_start
    isolde_start(session)
    from .place_ligand import place_ligand
    place_ligand(session, *args, **kwargs)

def add_water(session, *args, **kwargs):
    from chimerax.isolde.cmd import isolde_start
    isolde_start(session)
    from .place_ligand import place_water
    place_water(session, *args, **kwargs)



def register_building_commands(logger):
    register_isolde_replace(logger)
    register_isolde_add(logger)

def register_isolde_replace(logger):
    from chimerax.core.commands import (
        register, CmdDesc, StringArg
    )
    from chimerax.atomic import ResiduesArg
    desc = CmdDesc(
        required=[
            ('residue', ResiduesArg),
            ('new_residue_name', StringArg)
        ],
        synopsis='(EXPERIMENTAL) Replace an existing ligand with a related one, keeping as much of the original as possible'
    )
    register('isolde replace ligand', desc, replace_residue, logger=logger)

def register_isolde_add(logger):
    from chimerax.core.commands import (
        register, CmdDesc,
        AtomSpecArg, ModelArg,
        FloatArg, IntArg, BoolArg, StringArg, NoArg,
        ListOf, EnumOf, RepeatOf
    )
    from chimerax.atomic import AtomsArg, ResiduesArg, StructureArg

    def register_isolde_add_water():
        desc = CmdDesc(
            optional = [
                ('model', StructureArg),
            ],
            keyword = [
                ('position', ListOf(FloatArg,3,3)),
                ('bfactor', FloatArg),
                ('chain', StringArg),
                ('distance_cutoff', FloatArg),
                ('sim_settle', BoolArg)
            ],
            synopsis = 'Add a water molecule'
        )
        register('isolde add water', desc, add_water, logger=logger)
    register_isolde_add_water()

    def register_isolde_add_ligand():
        desc = CmdDesc(
            required = [
                ('ligand_id', StringArg),
            ],
            optional = [
                ('model', StructureArg),
            ],
            keyword = [
                ('position', ListOf(FloatArg,3,3)),
                ('bfactor', FloatArg),
                ('chain', StringArg),
                ('distance_cutoff', FloatArg),
                ('sim_settle', BoolArg),
                ('use_md_template', BoolArg),
                ('md_template_name', StringArg)
            ],
            synopsis = 'Add a ligand from the Chemical Components Dictionary'
        )
        register('isolde add ligand', desc, add_ligand, logger=logger)
    register_isolde_add_ligand()
