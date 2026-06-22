# @Author: Tristan Croll <tic20>
# @Date:   16-Apr-2020
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 05-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

from chimerax.isolde.cmd import block_if_sim_running

def replace_residue(session, residue, with_residue=None, apply_md_template=False):
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
    # The replacement template is a ChemComp identifier (a CCD code or a
    # registered id), resolved from the local store. Using an extant residue as a
    # one-off template was retired -- register it first with "isolde register
    # ligand", then replace by identifier.
    if not isinstance(with_residue, str):
        raise UserError(
            'isolde replace ligand requires a ChemComp identifier (string) for '
            '"withResidue". To use a built residue as a template, register it '
            'first with "isolde register ligand".')
    from .place_ligand import _ccd_template_residue, _ligand_residue_name
    template = _ccd_template_residue(session, with_residue)
    try:
        from ..template_utils import (
            fix_residue_from_template,
            fix_residue_to_match_md_template
        )
        corrected, unfixable = fix_residue_from_template(residue, template,
            rename_residue=False, match_by='element')
        # Name the rebuilt residue like a freshly-placed ligand (CCD code direct,
        # else LIGnn + isolde_chemcomp_id), so a long/registered id never becomes
        # the residue name.
        resname, chemcomp_id = _ligand_residue_name(residue.structure, with_residue)
        residue.name = resname
        if chemcomp_id is not None:
            residue.isolde_chemcomp_id = chemcomp_id
        if apply_md_template:
            ff = session.isolde.forcefield_mgr[session.isolde.sim_params.forcefield]
            ligand_db = session.isolde.forcefield_mgr.ligand_db(session.isolde.sim_params.forcefield)
            from chimerax.isolde.openmm.openmm_interface import find_residue_templates
            tdict = find_residue_templates(Residues([residue]), ff, ligand_db=ligand_db, logger=session.logger)
            md_template_name = tdict.get(0)
            if md_template_name is not None:
                fix_residue_to_match_md_template(session, residue, ff._templates[md_template_name], cif_template=template)
    finally:
        template.structure.delete()
    # Only speak up about chirality when something actually happened: report the
    # in-place corrections, and highlight only centres that could NOT be fixed
    # automatically (which genuinely need the user's attention).
    if corrected:
        session.logger.info('Corrected {} inverted chiral centre(s) ({}) in {}; '
            'run a quick energy minimisation to relax them.'.format(
                len(corrected), ','.join(corrected), residue.name))
    if unfixable:
        from chimerax.atomic import Atoms, Atom
        bad = Atoms([a for a in (residue.find_atom(n) for n in unfixable)
                     if a is not None])
        warn_str = ('Rebuilt ligand {} has chiral centre(s) at {} (highlighted) '
            'that could not be corrected automatically (no rotatable substituent). '
            'Check these carefully; if wrong, delete with "del #{}/{}:{}{}" and '
            'replace with "isolde add ligand {}".').format(
                residue.name, ','.join(unfixable),
                residue.structure.id_string, residue.chain_id, residue.number,
                residue.insertion_code, residue.name)
        session.selection.clear()
        bad.selected = True
        bad.draw_modes = Atom.BALL_STYLE
        from chimerax.isolde.view import focus_on_selection
        focus_on_selection(session, bad)
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



def add_registered_ligand_from_residue(
    session, residue, name=None, collection="user", sessionOnly=False
):
    block_if_sim_running(session)
    from chimerax.core.errors import UserError
    if len(residue) != 1:
        raise UserError('Select a single residue to register as a ligand template.')
    from .place_ligand import register_ligand
    cid = register_ligand(
        session,
        residue=residue[0],
        name=name,
        collection=collection,
        session_only=sessionOnly
    )
    session.logger.info(
        'Registered ligand template "%s". Place copies with "isolde add ligand %s".' %
        (cid, cid)
    )


def add_registered_ligand_from_smiles(
    session, smiles, name=None, collection="user", sessionOnly=False
):
    block_if_sim_running(session)
    from .place_ligand import register_ligand
    cid = register_ligand(
        session, smiles=smiles, name=name, collection=collection, session_only=sessionOnly
    )
    session.logger.info(
        'Registered ligand template "%s" from SMILES. Place copies with '
        '"isolde add ligand %s".' % (cid, cid)
    )


def register_building_commands(logger):
    register_isolde_replace(logger)
    register_isolde_add(logger)
    register_isolde_register_ligand(logger)


def register_isolde_register_ligand(logger):
    from chimerax.core.commands import register, CmdDesc, StringArg, BoolArg
    from chimerax.atomic import ResiduesArg
    res_desc = CmdDesc(
        required=[('residue', ResiduesArg)],
        keyword=[('name', StringArg), ('collection', StringArg), ('sessionOnly', BoolArg)],
        synopsis=(
            'Register a single residue as a reusable ligand template in the '
            'local ChemComp store (persistent by default; sessionOnly true '
            'keeps it only in this .cxs session).'
        )
    )
    register(
        'isolde register ligand', res_desc, add_registered_ligand_from_residue, logger=logger
    )
    smi_desc = CmdDesc(
        required=[('smiles', StringArg)],
        keyword=[('name', StringArg), ('collection', StringArg), ('sessionOnly', BoolArg)],
        synopsis=(
            'Register a ligand template from a SMILES string (a 3D '
            'conformer is generated) into the local ChemComp store. '
            'Requires name <identifier>.'
        )
    )
    register(
        'isolde register ligand smiles',
        smi_desc,
        add_registered_ligand_from_smiles,
        logger=logger
    )


def register_isolde_replace(logger):
    from chimerax.core.commands import (register, CmdDesc, StringArg, BoolArg)
    from chimerax.atomic import ResiduesArg
    desc = CmdDesc(
        required=[
            ('residue', ResiduesArg),
        ],
        keyword=[('with_residue', StringArg), ('apply_md_template', BoolArg)],
        synopsis=
        '(EXPERIMENTAL) Replace an existing ligand with a related one (by ChemComp identifier), keeping as much of the original as possible'
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
                # A ChemComp identifier: a CCD code, or an id registered with
                # "isolde register ligand" (resolved from the local ChemComp
                # store, with network fallback for CCD codes). Using an extant
                # residue as a one-off template was retired -- register it first.
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
            synopsis = 'Add a ligand by ChemComp identifier (CCD code or a registered id)'
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
        register('isolde modify his', desc, mod_his_command, logger=logger)
    register_isolde_mod_his()
