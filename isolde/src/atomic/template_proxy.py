# @Author: Tristan Croll
# @Date:   20-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 20-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

'''
A uniform "coordinate template" interface over either a ChimeraX ``TmplResidue``
(a CCD component) or a plain ``Residue``.

ISOLDE's residue fix/rebuild code historically required a ``TmplResidue`` as the
source of correct geometry, which means it can only repair residues that match a
CCD entry -- ChimeraX exposes no Python API to mint a *custom* ``TmplResidue``.
A user with a correctly-built but non-CCD ligand (e.g. loaded from a file) had no
way to use it as the repair template.

``TemplateProxy`` wraps either source behind one small interface (``name``,
``atoms``, ``find_atom``) and adds the RDKit chemistry the template needs for
chirality-/bond-order-aware matching via :func:`to_rdkit`. Because a plain
``Residue`` already provides the atom-access API the fix/rebuild code uses
(``atoms``, ``find_atom``, per-atom ``name``/``element``/``coord``/``neighbors``/
``bonds``), the proxy is deliberately thin -- its real value is the RDKit bridge.
'''

from . import rdkit_bridge as rb


class TemplateProxy:
    '''Wrap a ``TmplResidue`` or a ``Residue`` as a coordinate template.

    Args:
        session: the ChimeraX session.
        source: a :class:`chimerax.atomic.cytmpl.TmplResidue`, a
            :class:`chimerax.atomic.Residue`, or an existing ``TemplateProxy``
            (returned unchanged-in-spirit by unwrapping).
    '''

    def __init__(self, session, source):
        if isinstance(source, TemplateProxy):
            source = source.source
        self.session = session
        self._source = source
        self.name = source.name

    @property
    def source(self):
        '''The wrapped ``TmplResidue`` or ``Residue``.'''
        return self._source

    @property
    def is_residue(self):
        '''True if the wrapped source is a live ``Residue`` (vs a ``TmplResidue``).'''
        from chimerax.atomic import Residue
        return isinstance(self._source, Residue)

    @property
    def atoms(self):
        return self._source.atoms

    def find_atom(self, name):
        return self._source.find_atom(name)

    def to_rdkit(self):
        '''Return ``(mol, atom_map)`` for the template, with stereochemistry from
        coordinates and a ``cxName`` property per atom; ``atom_map`` maps RDKit
        atom index -> the wrapped template atom. Returns ``(None, {})`` on
        failure.

        * ``Residue`` source -> perceive/convert via
          :func:`rdkit_bridge.residue_to_rdkit` (which already returns the map).
        * ``TmplResidue`` source -> build from the CCD connection table
          (:func:`rdkit_bridge.template_to_rdkit`) and map back by ``cxName``.
        '''
        if self.is_residue:
            return rb.residue_to_rdkit(self._source)
        # Prefer the explicit CCD connection table (exact bond orders/charges);
        # fall back to perceiving from the template's own ideal geometry.
        mol, _status = rb.template_to_rdkit(self.session, self.name)
        if mol is None:
            mol = rb.template_geometry_mol(self._source)
        if mol is None:
            return None, {}
        atom_map = {}
        for a in mol.GetAtoms():
            if a.HasProp(rb.NAME_PROP):
                ta = self.find_atom(a.GetProp(rb.NAME_PROP))
                if ta is not None:
                    atom_map[a.GetIdx()] = ta
        return mol, atom_map
