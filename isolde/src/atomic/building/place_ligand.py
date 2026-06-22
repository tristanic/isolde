# @Author: Tristan Croll <tic20>
# @Date:   11-Jun-2019
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 04-Jan-2021
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll



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
    from chimerax.atomic.struct_edit import add_atom
    e = Element.get_element(element_name)
    r = model.new_residue(residue_name, chain_id, residue_number)
    add_atom(atom_name, e, r, coord, bfactor=bfactor)
    return r

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
    return place_ligand(session, 'HOH', model=model, position=position, bfactor=bfactor,
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
    from chimerax.geometry import find_closest_points
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
            bfactor = na.residue.atoms.bfactors.mean() + 5
    # The ligand may be a CCD id (built from the local ChemComp store) or a
    # user-supplied residue used as a custom geometry template. Resolve both to a
    # single coordinate-template residue; delete the throwaway structure (CCD path
    # only) once placement and any MD-template correction are done.
    tmpl, owns_template = _resolve_ligand_template(session, ligand_id)
    try:
        ligand_name = tmpl.name
        # Short, name-legal identifiers (e.g. CCD codes like NAG) become the
        # residue name directly; longer/registered identifiers get a short
        # generated LIGnn name, with the real identifier recorded on the residue
        # (isolde_chemcomp_id) so chirality/rebuild can re-find the chemistry.
        resname, chemcomp_id = _ligand_residue_name(model, ligand_name)
        r = new_residue_from_template(model, tmpl, chain, position,
                                      b_factor=bfactor, residue_name=resname)
        if chemcomp_id is not None:
            r.isolde_chemcomp_id = chemcomp_id
        if use_md_template and len(r.atoms) > 3:
            ff = session.isolde.forcefield_mgr[session.isolde.sim_params.forcefield]
            if md_template_name is None:
                from chimerax.isolde.openmm.amberff.template_utils import ccd_to_known_template
                md_template_name = ccd_to_known_template.get(ligand_name, None)
            if md_template_name is None:
                ligand_db = session.isolde.forcefield_mgr.ligand_db(
                    session.isolde.sim_params.forcefield
                )
                from chimerax.isolde.openmm.openmm_interface import find_residue_templates
                from chimerax.atomic import Residues
                tdict = find_residue_templates(
                    Residues([r]), ff, ligand_db=ligand_db, logger=session.logger
                )
                md_template_name = tdict.get(0)
            md_template = None
            if md_template_name is not None:
                md_template = ff._templates.get(md_template_name, None)

            if md_template is None:
                session.logger.warning(
                    'place_ligand() was called with use_md_template=True, '
                    'but no suitable template was found. This command has been ignored.'
                )
                return
            from ..template_utils import fix_residue_to_match_md_template
            fix_residue_to_match_md_template(session, r, md_template, cif_template=tmpl)
    finally:
        if owns_template:
            tmpl.structure.delete()


    matoms.selected=False
    r.atoms.selected=True
    if sim_settle:
        if not hasattr(session, 'isolde'):
            session.logger.warning('ISOLDE is not running. sim_settle argument ignored.')
        elif model != session.isolde.selected_model:
            session.logger.warning(
                "New ligand was not added to ISOLDE's "
                "selected model. sim_settle argument ignored."
            )
        else:
            # defer the simulation starting until after the new atoms have been
            # drawn, to make sure their styling "sticks"
            def do_run(*_, session=session):
                from chimerax.core.commands import run
                from chimerax.atomic import selected_residues, concise_residue_spec
                run(
                    session,
                    f'isolde sim start {concise_residue_spec(session, selected_residues(session))}'
                )
                from chimerax.core.triggerset import DEREGISTER
                return DEREGISTER

            session.triggers.add_handler('frame drawn', do_run)
    return r


def register_ligand(
    session, *, residue=None, smiles=None, name=None, collection="user", session_only=False
):
    """Register a ligand template -- from a modelled `residue` or a `smiles`
    string -- into the local ChemComp store under identifier `name`, so it can
    later be placed with ``isolde add ligand <name>``. Returns the stored id.

    Persistent by default (a user collection in its own ``<collection>.sqlite3``,
    safe from ``chemcomp update``); ``session_only=True`` stores it in the
    session collection that travels with the ``.cxs`` session file. ISOLDE owns
    the chemistry (records via :mod:`rdkit_bridge`); ChemComp only stores.
    """
    from chimerax.core.errors import UserError
    from ..rdkit_bridge import records_from_residue, records_from_smiles
    try:
        from chimerax import chemcomp
    except ImportError:
        raise UserError('The ChimeraX-ChemComp bundle is required to register ligands.')
    if residue is not None:
        rec = records_from_residue(residue, comp_id=(name or residue.name))
    elif smiles is not None:
        if not name:
            raise UserError('A name (identifier) is required when registering from SMILES.')
        rec = records_from_smiles(smiles, comp_id=name)
        if rec is None:
            raise UserError('Could not parse or embed SMILES: %r' % smiles)
    else:
        raise UserError('Provide a residue or a SMILES string to register.')
    return chemcomp.register_record(
        session, rec, collection=collection, persistent=not session_only
    )


def _ligand_residue_name(model, identifier):
    """Choose the in-model residue name for a ligand identified by `identifier`,
    returning ``(residue_name, chemcomp_id_to_record_or_None)``.

    * If a residue mapped to this identifier already exists in the model, reuse
      its name (so the same compound placed twice is named consistently).
    * A short, name-legal identifier (e.g. a CCD code) is used as the residue
      name directly -- no indirection needed.
    * Otherwise a fresh ``LIGnn`` name is minted and the real identifier is
      returned to be stored on the new residue's ``isolde_chemcomp_id``.
    """
    for r in model.residues:
        if getattr(r, 'isolde_chemcomp_id', None) == identifier:
            return r.name, (None if r.name == identifier else identifier)
        if r.name == identifier:
            return identifier, None
    if len(identifier) <= 5 and identifier.isalnum():
        return identifier, None
    existing = set(model.residues.names)
    i = 1
    while True:
        nm = 'LIG%02d' % i
        if nm not in existing:
            return nm, identifier
        i += 1


def _ccd_template_residue(session, ccd_id):
    '''Build a throwaway one-residue :class:`AtomicStructure` for a CCD component
    from ISOLDE's local ChemComp store, and return its :class:`Residue`.

    Uses :func:`rdkit_bridge.ccd_records`, which consults the local store first
    (instant, offline) and only falls back to the per-residue network fetch if the
    store misses -- so this is the store-backed replacement for
    ``mmcif.find_template_residue`` as a *geometry* template. The structure is not
    added to the session; the caller is responsible for deleting it
    (``residue.structure.delete()``).
    '''
    import numpy
    from chimerax.core.errors import UserError
    from chimerax.atomic import AtomicStructure, Element
    from chimerax.atomic.struct_edit import add_atom, add_bond
    from chimerax.isolde.atomic.rdkit_bridge import ccd_records
    recs = ccd_records(session, ccd_id)
    if recs is None:
        raise UserError(
            'Could not find a definition for "{}" in the local ChemComp store or '
            'the Chemical Components Dictionary. For a novel compound, register it '
            'first with "isolde register ligand".'.format(ccd_id)
        )
    atoms, bonds, coords = recs
    s = AtomicStructure(session, name='{} template'.format(ccd_id), auto_style=False)
    r = s.new_residue(ccd_id, 'A', 1)
    amap = {}
    for (aid, sym, _charge, _arom) in atoms:
        xyz = coords.get(aid)
        if xyz is None:
            # No ideal/model coordinate for this atom -- can't place it.
            continue
        e = Element.get_element((sym or 'C').strip().capitalize())
        amap[aid] = add_atom(aid, e, r, numpy.array(xyz, dtype=numpy.float64))
    for (a1, a2, _order, _arom) in bonds:
        x1, x2 = amap.get(a1), amap.get(a2)
        if x1 is not None and x2 is not None and x2 not in x1.neighbors:
            add_bond(x1, x2)
    return r


def _resolve_ligand_template(session, identifier):
    '''Build a coordinate-template residue for `identifier` (a CCD code or a
    registered ChemComp identifier) from the local ChemComp store, returning
    ``(template_residue, owns_structure=True)`` -- the caller deletes the
    throwaway structure.

    Using an extant residue as a one-off template was retired in favour of the
    register/apply split: register the residue once with ``isolde register
    ligand``, then place copies by identifier.
    '''
    from chimerax.core.errors import UserError
    if not isinstance(identifier, str):
        raise UserError(
            'isolde add ligand expects a ChemComp identifier (a string). To use a '
            'built residue as a template, register it first with "isolde register '
            'ligand".'
        )
    return _ccd_template_residue(session, identifier), True


def new_residue_from_template(
    model,
    template,
    chain_id,
    center,
    residue_number=None,
    insert_code=' ',
    b_factor=50,
    precedes=None,
    residue_name=None
):
    '''
    Create a new residue based on a template, and add it to the model. The new
    residue is named `residue_name` if given, else after the template.
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
    r = model.new_residue(residue_name or template.name, chain_id, residue_number,
        insert=insert_code, precedes=precedes)
    from chimerax.atomic.struct_edit import add_bond, add_atom
    for i, ta in enumerate(template.atoms):
        a = tatom_to_atom[ta] = add_atom(ta.name, ta.element, r, t_coords[i], bfactor=b_factor)
        for tn in ta.neighbors:
            n = tatom_to_atom.get(tn, None)
            if n is not None:
                add_bond(a, n)
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
