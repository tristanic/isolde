# @Author: Tristan Croll <tic20>
# @Date:   20-Jun-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 03-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

from chimerax.atomic import selected_residues, Residue
from chimerax.core.errors import UserError

tooltip = ('Add an OXT atom to a C-terminal amino acid residue.')

def add_oxt(session, residue):
    catom = residue.find_atom('C')
    resname = residue.name
    if catom is None:
        session.logger.warning('Residue {} {}{} has no C atom!'.format(residue.name, residue.chain_id, residue.number))
        return
    color = catom.color
    for n in catom.neighbors:
        if n.name=='OXT':
            session.logger.warning('Chain {} already has a C-terminal OXT. Skipping.')
            return
        if n.residue != residue:
            raise UserError('Residue {} {}{} is not a C-terminal residue!'.format(residue.name, residue.chain_id, residue.number))
    from chimerax.build_structure import modify_atom
    from chimerax.atomic import Element
    atoms = modify_atom(catom, catom.element, 3, res_name=residue.name)
    for a in atoms:
        if a.element.name=='H':
            break
    modify_atom(a, Element.get_element('O'), 1, name='OXT', res_name=residue.name)
    catom.color = color
    session.logger.info('Added a C-terminal OXT to chain {}'.format(residue.chain_id))

def run_script(session):
    sel = selected_residues(session)
    if not len(sel):
        raise UserError('No atoms selected!')
    sel = sel[sel.polymer_types==Residue.PT_AMINO]
    if not len(sel):
        raise UserError('No protein selected!')
    us = sel.unique_structures
    if len(us) != 1:
        raise UserError('Selection must be from a single model!')
    m = sel.unique_structures
    cids = sel.unique_chain_ids
    for cid in cids:
        chain = m.chains[m.chains.chain_ids==cid][0]
        last_res = chain.residues[-1]
        add_oxt(session, last_res)
