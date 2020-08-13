# @Author: Tristan Croll <tic20>
# @Date:   11-Jun-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 10-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

from . import set_new_atom_style


def _get_metals():
    '''
    Returns a list of metals sorted by mass
    '''
    from chimerax.atomic import Element
    metals = [n for n in Element.names if Element.get_element(n).is_metal]
    sorted_metals = sorted(metals, key=lambda m: Element.get_element(m).mass)
    return [Element.get_element(m) for m in sorted_metals]

_metals = _get_metals()

def place_metal_at_coord(model, chain_id, residue_number, residue_name, atom_name, coord, element_name=None, bfactor=20):
    '''
    Create a new residue encompassing a single metal ion in the current model,
    and place it at the given coordinate.
    '''
    if element_name is None:
        element_name = atom_name.title()
    from chimerax.atomic import Element
    e = Element.get_element(element_name)
    r = model.new_residue(residue_name, chain_id, residue_number)
    a = model.new_atom(atom_name, e)
    a.coord = coord
    a.bfactor=bfactor
    r.add_atom(a)

def place_water(session, model=None, position=None, bfactor=None, chain=None,
        distance_cutoff=3.0, sim_settle=True):
    '''
    Place a water molecule at the given position, or the current centre of
    rotation if no :attr:`position` is specified. If the :attr:`bfactor` or
    :attr:`chain` arguments are not specified, the water will be assigned the
    missing properties based on the nearest residue within
    :attr:`distance_cutoff` of :attr:`position`, with the assigned B-factor
    being (average B-factor of nearest residue)+5. If :attr:`sim_settle` is
    True, a small local simulation will automatically start to settle the
    water into position.
    '''
    place_ligand(session, 'HOH', model=model, position=position, bfactor=bfactor,
        chain=chain, distance_cutoff=distance_cutoff, sim_settle=sim_settle)

def place_ligand(session, ligand_id, model=None, position=None, bfactor=None, chain=None,
        distance_cutoff=8.0, sim_settle=False, use_md_template=True, md_template_name=None):
    '''
    Place a ligand at the given position, or the current centre of
    rotation if no :attr:`position` is specified. If the :attr:`bfactor` or
    :attr:`chain` arguments are not specified, the water will be assigned the
    missing properties based on the nearest residue within
    :attr:`distance_cutoff` of :attr:`position`, with the assigned B-factor
    being (average B-factor of nearest residue)+5. If :attr:`sim_settle` is
    True, a small local simulation will automatically start to settle the
    ligand into position.

    For residues of more than three atoms, if you have loaded an MD template for
    the residue (with matching residue name), any mis-matches between atoms in
    the CIF and MD templates will be automatically corrected as long as the
    changes only involve removing atoms from the CIF template and/or adding
    atoms that are directly connected to a single atom in the CIF template. The
    most common scenario where this arises is where the protonation state of the
    CIF template is different from that of the MD template.

    Note that at this stage no attempt is made at a preliminary fit to the
    density, or to avoid clashes: the ligand is simply placed at the given
    position with coordinates defined by the template in the Chemical
    Components Dictionary. Except for small, rigid ligands the use of
    :attr:`sim_settle`=True is inadvisable. If the ligand as placed clashes
    badly with the existing model, the best approach is to run the command
    `isolde ignore ~selAtoms` then start a simulation (which will involve *only*
    the new residue) to pull it into position (you may want to use pins where
    necessary). Once satisfied that there are no severe clashes, stop the
    simulation and run `isolde ~ignore` to reinstate all atoms for simulation
    purposes, and continue with your model building.
    '''
    from chimerax.core.geometry import find_closest_points
    from chimerax.atomic import mmcif
    from chimerax.core.errors import UserError
    if hasattr(session, 'isolde') and session.isolde.simulation_running:
        raise UserError('Cannot add atoms when a simulation is running!')
    if model is None:
        if hasattr(session, 'isolde') and session.isolde.selected_model is not None:
            model = session.isolde.selected_model
        else:
            raise UserError('If model is not specified, ISOLDE must be started '
                'with a model loaded')
    if position is None:
        position = session.view.center_of_rotation
    matoms = model.atoms
    if bfactor is None or chain is None:
        _,_,i = find_closest_points([position], matoms.coords, distance_cutoff)
        if not len(i):
            err_str = ('No existing atoms found within the given distance cutoff '
                'of the target position. You may repeat with a larger cutoff or '
                'explicitly specify the B-factor and chain ID')
            if ligand_id == 'HOH':
                err_str += (', but keep in mind that placing waters outside of '
                    'H-bonding distance to the model is generally inadvisable.')
            else:
                err_str += '.'
            raise UserError(err_str)
        na = matoms[i[0]]
        if chain is None:
            chain = na.residue.chain_id
        if bfactor is None:
            bfactor = na.residue.atoms.bfactors.mean()+5
    tmpl = mmcif.find_template_residue(session, ligand_id)
    r = new_residue_from_template(model, tmpl, chain, position, b_factor=bfactor)
    if use_md_template and len(r.atoms) > 3:
        ff = session.isolde.forcefield_mgr[session.isolde.sim_params.forcefield]
        if md_template_name is None:
            from chimerax.isolde.openmm.amberff.template_utils import ccd_to_known_template
            md_template_name = ccd_to_known_template.get(ligand_id, None)
        if md_template_name is None:
            ligand_db = session.isolde.forcefield_mgr.ligand_db(session.isolde.sim_params.forcefield)
            from chimerax.isolde.openmm.openmm_interface import find_residue_templates
            from chimerax.atomic import Residues
            tdict = find_residue_templates(Residues([r]), ff, ligand_db=ligand_db, logger=session.logger)
            md_template_name = tdict.get(0)
        md_template = None
        if md_template_name is not None:
            md_template = ff._templates.get(md_template_name, None)

        if md_template is None:
            session.logger.warning('place_ligand() was called with use_md_template=True, '
                'but no suitable template was found. This command has been ignored.')
            return
        from ..template_utils import fix_residue_to_match_md_template
        fix_residue_to_match_md_template(session, r, md_template, cif_template=tmpl)


    matoms.selected=False
    r.atoms.selected=True
    if sim_settle:
        if not hasattr(session, 'isolde'):
            session.logger.warning('ISOLDE is not running. sim_settle argument ignored.')
        elif model != session.isolde.selected_model:
            session.logger.warning("New ligand was not added to ISOLDE's "
                "selected model. sim_settle argument ignored.")
        else:
            from chimerax.core.commands import run
            run(session, 'isolde sim start sel')

def new_residue_from_template(model, template, chain_id, center,
        residue_number=None, insert_code=' ', b_factor=50, precedes=None):
    '''
    Create a new residue based on a template, and add it to the model.
    '''
    if residue_number is None:
        if chain_id in model.residues.chain_ids:
            residue_number = suggest_new_residue_number_for_ligand(model, chain_id)
        else:
            residue_number = 0
    import numpy
    from chimerax.atomic import Atom
    t_coords = numpy.array([a.coord for a in template.atoms])
    t_center = t_coords.mean(axis=0)
    t_coords += numpy.array(center) - t_center
    tatom_to_atom = {}
    r = model.new_residue(template.name, chain_id, residue_number,
        insert=insert_code, precedes=precedes)
    for i, ta in enumerate(template.atoms):
        a = tatom_to_atom[ta] = model.new_atom(ta.name, ta.element)
        a.coord = t_coords[i]
        a.bfactor = b_factor
        r.add_atom(a)
        for tn in ta.neighbors:
            n = tatom_to_atom.get(tn, None)
            if n is not None:
                model.new_bond(a, n)
    set_new_atom_style(model.session, r.atoms)
    return r



def find_nearest_chain(model, coord):
    '''
    Find the chain coming nearest to the given coordinate.
    '''
    pass

def suggest_new_residue_number_for_ligand(model, chain_id):
    '''
    Suggest a suitable residue number for a new ligand, based on what is already
    in the chain.
    '''
    from chimerax.atomic import Residue
    residues = model.residues[model.residues.chain_ids == chain_id]
    ligand_residues = residues[residues.polymer_types==Residue.PT_NONE]
    if ligand_residues:
        return max(ligand_residues.numbers)+1

    last_polymeric_residue_num = max(residues.numbers)
    # Start ligands on a nice round multiple of 1000
    ligand_num = round(last_polymeric_residue_num+1000, -3)
    if ligand_num > 9999:
        raise TypeError('The PDB format does not support residue numbers greater '
            'than 9999. Consider adding your ligand to a different chain.')
    return ligand_num
