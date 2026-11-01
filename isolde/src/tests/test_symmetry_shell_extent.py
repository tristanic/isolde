# @Author: Tristan Croll
# @Date:   29-Oct-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 29-Oct-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Regression test for how far a crystallographic symmetry shell and its map
coverage extend.

  1. Fragment completion (``symmetry_sim._shell_fragments``): a protein shell is
     unchanged (whole residues), but a small-molecule (corecif) model -- whose
     whole asymmetric unit is one residue -- is completed per molecule. The
     special-position check is per molecule too: in C2/c COD 1100908 the Cu sits
     on a 2-fold, and the old per-residue check pruned the whole residue's image,
     dropping the two ligand images that complete the Cu coordination sphere.
  2. Map coverage (``SimManager._cover_maps_including_symmetry``): only the
     region of interest, its immediate crystal contacts, and the one image each
     promoted parent's MDFF term samples are covered -- not every shell copy,
     which made the live map box the bounding box of the whole shell.
  3. No duplicate images (``symmetry_sim._duplicate_fragment_images``): a molecule
     or atom with site symmetry ``S`` makes ``g.M`` and ``g.S.M`` the same neighbour,
     so the shell used to build it twice, atom on top of atom. For a molecule held
     whole in the model (``clipper fragments complete``) the coincident twins got a
     nonzero copy<->copy weight and the sim exploded at once (COD 2201140).

Uses Clipper's bundled test crystals and lightweight fakes for the map manager,
so it needs no running simulation or GUI. Run inside ChimeraX:

    run_chimerax.bat --nogui --exit --script src/tests/test_symmetry_shell_extent.py

Prints PASS/FAIL and exits non-zero on failure.
'''
import os
import numpy
from types import SimpleNamespace


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


def _clipper_test_file(name):
    import chimerax.clipper
    return os.path.join(os.path.dirname(chimerax.clipper.__file__), 'tests', name)


def _copies_by_op(shell):
    from chimerax.atomic import Atoms
    out = {}
    for c in shell.copies:
        out.setdefault(c.symop_index, []).append(c.parent_atom)
    return {k: Atoms(v) for k, v in out.items()}


def _ghosts(shell, op, atoms):
    m = shell.symmats[op]
    return atoms.coords @ m[:, :3].T + m[:, 3]


def test_protein_shell_whole_residues(session):
    from chimerax.core.commands import run
    from chimerax.isolde.openmm.symmetry_sim import build_symmetry_shell
    p = run(session, 'open ' + _clipper_test_file('3io0.pdb'), log=False)[0]
    try:
        shell = build_symmetry_shell(p, p.residues[:40].atoms, 8.0)
        if shell is None or not len(shell.copies):
            _fail('no symmetry copies for the 3io0 test selection')
        for op, atoms in _copies_by_op(shell).items():
            res = atoms.unique_residues
            if len(res.atoms) != len(atoms):
                _fail('3io0: operator %d copies partial residues' % op)
        print('PASS: protein shell is whole residues (%d copies)' % len(shell.copies))
    finally:
        session.models.close([p])


def test_small_molecule_shell_per_molecule(session):
    from chimerax.clipper.io.small_molecule import open_small_molecule_cif
    from chimerax.isolde.openmm.symmetry_sim import build_symmetry_shell
    m = open_small_molecule_cif(session, _clipper_test_file('cod_1100908.cif'))
    try:
        if m.num_residues != 1 or len(m.molecules) < 2:
            _fail('expected one corecif residue holding several molecules')
        shell = build_symmetry_shell(m, m.atoms, 12.0)
        whole_asu = len(m.atoms) * len(_copies_by_op(shell))
        if len(shell.copies) >= whole_asu:
            _fail('copies are not per molecule (%d copies, %d if whole ASU)'
                  % (len(shell.copies), whole_asu))
        mol_of = numpy.empty(len(m.atoms), dtype=int)
        for i, mol in enumerate(m.molecules):
            mol_of[m.atoms.indices(mol)] = i
        sizes = numpy.bincount(mol_of)
        cu = m.atoms[m.atoms.element_names == 'Cu']
        cu_xyz = cu.coords[0]
        coordinating = 0
        for op, atoms in _copies_by_op(shell).items():
            mols = mol_of[m.atoms.indices(atoms)]
            for k in numpy.unique(mols):
                if (mols == k).sum() != sizes[k]:
                    _fail('operator %d copies part of a molecule' % op)
            g = _ghosts(shell, op, atoms)
            d = numpy.linalg.norm(g - cu_xyz, axis=1)
            if (atoms.element_names == 'Cu').any() and d[atoms.element_names == 'Cu'].min() < 0.5:
                _fail('Cu image coincident with the real Cu under operator %d was kept' % op)
            heavy = atoms.element_names != 'H'
            coordinating += int(((d < 2.1) & heavy & (atoms.element_names != 'Cu')).sum())
        if coordinating < 2:
            _fail('symmetry-mate ligands coordinating the real Cu are missing (%d donor atoms)'
                  % coordinating)
        print('PASS: small-molecule shell is per molecule (%d copies vs %d whole-ASU); '
              '%d symmetry-mate donor atoms coordinate the real Cu'
              % (len(shell.copies), whole_asu, coordinating))
    finally:
        session.models.close([m])


class _FakeMapMgr:
    def cover_atoms(self, atoms, transforms=None, transform_indices=None,
            mask_radius=3, extra_padding=12):
        self.atoms, self.transforms, self.tidx = atoms, transforms, transform_indices


def _covered(model, roi, shell, mask_radius=4.0):
    '''Run the real coverage method on fakes; return (covered copy keys, fake).'''
    from chimerax.isolde.openmm.openmm_interface import SimManager
    mgr = _FakeMapMgr()
    fake_self = SimpleNamespace(_symmetry_shell=shell, model=model,
        sim_construct=SimpleNamespace(mobile_atoms=roi))
    SimManager._cover_maps_including_symmetry(fake_self,
        SimpleNamespace(map_mgr=mgr), SimpleNamespace(map_mask_radius=mask_radius), roi)
    n_roi = int((roi.element_names != 'H').sum())
    mats = mgr.transforms.array()
    keys = set()
    for a, ti in zip(mgr.atoms[n_roi:], mgr.tidx[n_roi:]):
        op = next(k for k in range(len(shell.symmats))
                  if numpy.allclose(shell.symmats[k], mats[ti], atol=1e-6))
        keys.add((int(model.atoms.index(a)), op))
    return keys, mgr


def test_coverage_whole_model(session):
    from chimerax.clipper.io.small_molecule import open_small_molecule_cif
    from chimerax.isolde.openmm.symmetry_sim import build_symmetry_shell
    from chimerax.geometry import find_close_points
    m = open_small_molecule_cif(session, _clipper_test_file('cod_1100908.cif'))
    try:
        shell = build_symmetry_shell(m, m.atoms, 12.0)
        keys, mgr = _covered(m, m.atoms, shell)
        if not keys or len(keys) >= len(shell.copies) / 2:
            _fail('whole-model coverage not trimmed (%d of %d copies covered)'
                  % (len(keys), len(shell.copies)))
        heavy = m.atoms[m.atoms.element_names != 'H'].coords
        for ai, op in keys:
            g = _ghosts(shell, op, m.atoms[ai:ai + 1])
            i1, _ = find_close_points(g, heavy, 4.0 + 1e-6)
            if not len(i1):
                _fail('covered a copy that is not a contact of the region of interest')
        print('PASS: whole-model coverage trimmed to %d contact copies (of %d)'
              % (len(keys), len(shell.copies)))
    finally:
        session.models.close([m])


def test_coverage_promoted_parents(session):
    from chimerax.core.commands import run
    from chimerax.isolde.openmm.symmetry_sim import build_symmetry_shell
    p = run(session, 'open ' + _clipper_test_file('3io0.pdb'), log=False)[0]
    try:
        roi = p.residues[:40].atoms
        shell = build_symmetry_shell(p, roi, 8.0)
        keys, mgr = _covered(p, roi, shell)
        centroid = roi.coords.mean(axis=0)
        images = {}
        for c in shell.copies:
            if roi.index(c.parent_atom) >= 0:
                continue
            m = shell.symmats[c.symop_index]
            g = m[:, :3] @ c.parent_atom.coord + m[:, 3]
            images.setdefault(int(p.atoms.index(c.parent_atom)), []).append(
                (c.symop_index, g))
        if not images:
            _fail('test selection has no promoted parents - not exercising anything')
        for ai, imgs in images.items():
            # The image _add_mdff_symmetry_terms gives this atom's map term.
            op = min(imgs, key=lambda im: numpy.sum((im[1] - centroid) ** 2))[0]
            if (ai, op) not in keys:
                _fail('promoted parent atom %d: its MDFF image (operator %d) is not '
                      'covered' % (ai, op))
        print('PASS: every promoted parent (%d atoms) has its MDFF image covered; '
              '%d of %d copies covered' % (len(images), len(keys), len(shell.copies)))
    finally:
        session.models.close([p])


def _coincidences(model, shell):
    '''(copy pairs closer than 0.05 A, copies within 0.5 A of a real atom).'''
    from scipy.spatial import cKDTree
    from chimerax.atomic import Atoms
    par = Atoms([c.parent_atom for c in shell.copies])
    ops = numpy.array([c.symop_index for c in shell.copies])
    mats = numpy.asarray(shell.symmats)[ops]
    g = numpy.einsum('nij,nj->ni', mats[:, :, :3], par.coords) + mats[:, :, 3]
    pairs = len(cKDTree(g).query_pairs(0.05))
    on_real = int(numpy.isfinite(cKDTree(model.atoms.coords).query(
        g, distance_upper_bound=0.5)[0]).sum())
    return pairs, on_real


def test_no_duplicate_images(session):
    from chimerax.core.commands import run
    from chimerax.clipper.io.small_molecule import open_small_molecule_cif
    from chimerax.clipper.cmd import split_fragments_cmd
    from chimerax.isolde.openmm.symmetry_sim import build_symmetry_shell
    cases = [
        # (file, fragments mode, what has site symmetry)
        ('cod_2021365.cif', 'complete', 'molecule completed across an inversion centre (P-1)'),
        ('cod_2237730.cif', 'rename', 'molecule lying in a mirror plane (P2_1/m)'),
        ('cod_1100908.cif', 'rename', 'Cu on a 2-fold (C2/c)'),
    ]
    for name, mode, what in cases:
        m = open_small_molecule_cif(session, _clipper_test_file(name))
        try:
            n0 = m.num_atoms
            split_fragments_cmd(session, [m], mode=mode)
            if mode == 'complete' and m.num_atoms <= n0:
                _fail('%s: fragments complete added no atoms - not exercising a whole '
                      'molecule on a special position' % name)
            shell = build_symmetry_shell(m, m.atoms, 12.0)
            pairs, on_real = _coincidences(m, shell)
            if pairs or on_real:
                _fail('%s (%s): %d coincident copy pairs, %d copies on real atoms'
                      % (name, what, pairs, on_real))
            print('PASS: %s, %s: %d copies, none coincident' % (name, what, len(shell.copies)))
        finally:
            session.models.close([m])
    p = run(session, 'open ' + _clipper_test_file('3io0.pdb'), log=False)[0]
    try:
        shell = build_symmetry_shell(p, p.atoms, 8.0)
        pairs, on_real = _coincidences(p, shell)
        if pairs or on_real:
            _fail('3io0: %d coincident copy pairs, %d copies on real atoms' % (pairs, on_real))
        print('PASS: 3io0 protein shell (%d copies) has no coincident images' % len(shell.copies))
    finally:
        session.models.close([p])


def run_all(session):
    test_protein_shell_whole_residues(session)
    test_small_molecule_shell_per_molecule(session)
    test_coverage_whole_model(session)
    test_coverage_promoted_parents(session)
    test_no_duplicate_images(session)
    print('ALL PASS')


if 'session' in globals():
    run_all(session)
