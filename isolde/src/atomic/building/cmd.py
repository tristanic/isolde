# @Author: Tristan Croll <tic20>
# @Date:   16-Apr-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 05-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

from chimerax.isolde.cmd import block_if_sim_running

def replace_residue(session, residue, new_residue_name):
    block_if_sim_running(session)
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
    from chimerax import mmcif
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
    block_if_sim_running(session)
    from chimerax.isolde.cmd import isolde_start
    isolde_start(session)
    from .place_ligand import place_ligand
    place_ligand(session, *args, **kwargs)

def add_water(session, *args, **kwargs):
    block_if_sim_running(session)
    from chimerax.isolde.cmd import isolde_start
    isolde_start(session)
    from .place_ligand import place_water
    place_water(session, *args, **kwargs)

def add_aa_cmd(session, resname, stem_residue=None, add_direction=None, structure=None, chain_id=None, number=None, add_b_factor=0, occupancy=1,
        approx_conformation='strand'):
    block_if_sim_running(session)
    from chimerax.core.errors import UserError
    if stem_residue is None and (structure is None or chain_id is None or number is None):
        raise UserError('If stem residue is not specified, structure, chain ID and number must be provided!')
    if stem_residue is None:
        residue = None
    else:
        if len(stem_residue) != 1:
            raise UserError('Please select a single residue!')
        residue = stem_residue[0]
        structure = residue.structure
        from chimerax.atomic import Residue
        if residue.polymer_type != Residue.PT_AMINO:
            raise UserError('Selection must be an amino acid from a protein chain!')
    if approx_conformation=='strand':
        phi=-139
        psi=135
    elif approx_conformation=='helix':
        phi=-64
        psi=-41
    else:
        raise UserError('approx_conformation must be one of "strand" or "helix"!')
    prev_res = next_res = None
    if residue is not None:
        aa_neighbors = []
        for n in residue.neighbors:
            if n.polymer_type == Residue.PT_AMINO:
                is_peptide_bond = False
                for b in residue.bonds_between(n):
                    if all([a.element.name in ('N','C') for a in b.atoms]):
                        is_peptide_bond = True
                if is_peptide_bond:
                    aa_neighbors.append(n)
        if len(aa_neighbors) > 1:
            raise UserError('Selection is not a terminal residue!')
        elif len(aa_neighbors) == 0: 
            if add_direction is None:
                raise UserError('Stem residue has no existing neighbors. Please specify a build direction')
            if add_direction == 'N':
                next_res = residue
            elif add_direction == 'C':
                prev_res = residue
            else:
                raise UserError('Stem residue must be one of "N" or "C"!')
        else: 
            n = aa_neighbors[0]
            if (n.number < residue.number or 
                (n.number == residue.number and n.insertion_code < residue.insertion_code)):
                prev_res = residue
            else:
                next_res = residue
    from .build_utils import add_amino_acid_residue
    add_amino_acid_residue(structure, resname.upper(), prev_res=prev_res, next_res=next_res,
        chain_id=chain_id, number=number, add_b_factor=add_b_factor, occupancy=occupancy, phi=phi, psi=psi)

def mod_his_command(session, residues, protonated_atom):
    block_if_sim_running(session)
    residues = residues[residues.names=='HIS']
    if not len(residues):
        from chimerax.core.errors import UserError
        raise UserError('No histidine residues specified!')
    from .build_utils import set_his_protonation_state
    for r in residues:
        set_his_protonation_state(r, protonated_atom)
        session.logger.info(f'Set protonation of HIS #{r.structure.id_string}/{r.chain_id}:{r.number}{r.insertion_code} to {protonated_atom}')



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
        ListOf, EnumOf, RepeatOf, Bounded
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

    def register_isolde_add_aa():
        desc = CmdDesc(
            required=[
                ('resname', StringArg),
            ],
            optional=[
                ('stem_residue', ResiduesArg),
            ],
            keyword=[
                ('add_direction', EnumOf(['C', 'N'])),
                ('structure', StructureArg),
                ('chain_id', StringArg),
                ('number', IntArg),
                ('add_b_factor', Bounded(FloatArg, min=0)),
                ('occupancy', Bounded(FloatArg, min=0, max=1)),
                ('approx_conformation', EnumOf(['strand','helix']))
            ],
            synopsis='Add an amino acid to the terminus of an existing chain, or as the start of a new fragment.'
        )
        register('isolde add aa', desc, add_aa_cmd, logger=logger)
    register_isolde_add_aa()
    
    def register_isolde_mod_his():
        desc = CmdDesc(
            synopsis='Set the protonation state of one or more histidine residues',
            required = [
                ('residues', ResiduesArg),
                ('protonated_atom', EnumOf(['ND', 'NE', 'both']))
            ]
        )
        register('isolde mod his', desc, mod_his_command, logger=logger)
    register_isolde_mod_his()
