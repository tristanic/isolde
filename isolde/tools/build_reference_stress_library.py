# @Author: Tristan Croll
# @Date:   01-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 01-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

'''
Dev-time builder for validation/reference_stress.ReferenceStressLibrary's baseline
cache (validation/data/reference_virial.npz). NOT part of the shipped bundle --
run this once, interactively, whenever the reference set or force field changes,
and check the resulting .npz into the repo.

Like isolde/src/tests/test_simulation.py, this is a manual, interactively-run tool,
not an automated test -- ISOLDE has no pytest/CI (see the project's CLAUDE.md).
`tools/` is a dev-only directory, not part of the installed chimerax.isolde
package (cf. tools/gen_rs_glyphs.py's own docstring) -- so it can't be imported
by dotted path. Run it from ChimeraX's Python shell (Tools -> General -> Shell)
with `session` already bound, loading the file directly by path:

    import importlib.util
    spec = importlib.util.spec_from_file_location(
        'build_reference_stress_library',
        '/path/to/isolde/isolde/tools/build_reference_stress_library.py')
    brsl = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(brsl)
    brsl.build(session)

This fetches each reference structure from the PDB, builds an ISOLDE simulation
over it (via SimHandler._prepare_sim(), the same call start_sim() itself makes
first -- see openmm_interface.py:1549 -- but skipping the rest of start_sim()'s
async, GUI-frame-triggered continuation logic, which is designed for live
interactive stepping rather than a single deterministic minimize-and-extract
pass), minimizes once at ISOLDE's own default tolerance, and accumulates a
per-(residue name, atom name) sample of the symmetrised virial tensor across all
structures. NOTE: every individual API call below was confirmed against the
current openmm_interface.py source (not copied from the possibly-stale
test_simulation.py, whose SimHandler/SimConstruct call omits an argument the
current signature requires) -- but the script as a whole has not been
execution-tested end to end, since that requires a live GUI ChimeraX+ISOLDE
session with real network fetches. Sanity-check its output before trusting it.
'''

import numpy as np

# Curated, well-refined reference structures spanning protein + nucleic acid.
# Each entry is a dict:
#   'id'         : PDB ID (required)
#   'assembly'   : (optional) mmCIF author-defined assembly id to expand+bake into
#                  a single structure before simulating -- needed when the AU is
#                  only *part* of the biological unit (e.g. a single DNA strand
#                  whose duplex partner is a symmetry mate). Without this the lone
#                  fragment minimises with no base-pairing partner -> unrealistic
#                  baseline.
#   'delete'     : (optional) extra atomspec fragment deleted after cleaning, e.g.
#                  a modified residue ISOLDE can't parameterise from standard
#                  templates. Deleting mid-chain leaves a break (harmless -- ISOLDE
#                  builds with ignoreExternalBonds=True), but see the note in
#                  _minimize_and_get_virial about the small terminal-contamination
#                  it introduces.
# Not claimed to be complete or uniformly sub-Angstrom -- extend via a
# resolution-filtered RCSB search. Glycan references are still future work.
DEFAULT_STRUCTURES = [
    {'id': '1EJG'},                    # crambin, 0.54 A -- protein
    {'id': '2VB1'},                    # hen egg-white lysozyme, ~0.65 A -- protein
    {'id': '2IXT'},                    # sphericase -- larger protein example
    {'id': '3GGK', 'assembly': '1'},   # B-DNA; AU is one strand, assembly 1 (x2)
                                       # builds the duplex
    {'id': '1Y26'},                    # adenine riboswitch aptamer -- RNA, 71 nt
                                       # single chain, unmodified A/C/G/U (bound
                                       # adenine ligand is stripped as non-polymer)
    {'id': '1MSY'},                    # loop E of 5S rRNA, 1.55 A -- clean RNA, 27 nt
    # 5MRX (RNA) was tried but doesn't parameterise headlessly: it combines a
    # modified nucleotide (5HM) with internal chain breaks, and each resulting free
    # 5'-OH segment end gets an HO5' hydrogen ISOLDE's AMBER RNA templates reject
    # (unmatched: U2647/U2650/A2654). Standard RNA termini are fine -- 1Y26/1MSY
    # above parameterise cleanly -- so this is specific to 5MRX's chemistry, not an
    # RNA-general limitation. Revisit only with proper terminal/chain-break handling.
]


def _normalize_spec(spec):
    return {'id': spec} if isinstance(spec, str) else dict(spec)


def _is_terminal_residue(residue):
    chain = residue.chain
    if chain is None:
        return True  # not part of a polymer chain at all -- treat as an edge case
    residues = chain.residues
    return residue is residues[0] or residues[-1] is residue


def _prepare_structure(session, spec):
    '''Fetch spec['id'] (optionally expanding an author-defined assembly), strip
    solvent/heteroatoms, and return the cleaned single AtomicStructure. Keeps
    standard protein/nucleic polymer only -- ligands/novel residues need their own
    parameterisation path and aren't part of this baseline.'''
    from chimerax.core.commands import run
    from chimerax.atomic import AtomicStructure
    pdb_id = spec['id']
    assembly = spec.get('assembly')

    au = run(session, f'open {pdb_id}')[0]
    if assembly is None:
        model = au
    else:
        # `sym ... copies true` makes one model copy per assembly operator, each
        # sharing the AU's local coords but carrying its own position transform.
        # `combine` then bakes those transforms into a single structure in a common
        # coordinate frame (remapping colliding chain IDs) -- which is what ISOLDE
        # needs, since SimConstruct reads atoms.coords directly and ignores per-model
        # position transforms.
        before = set(session.models.list())
        run(session, f'sym #{au.id_string} assembly {assembly} copies true')
        after = set(session.models.list())
        copies = [m for m in (after - before) if isinstance(m, AtomicStructure)]
        if not copies:
            raise RuntimeError(f'assembly {assembly!r} produced no atomic copies')
        from chimerax.atomic.cmd import combine_cmd
        model = combine_cmd(session, copies, close=True,
                            name=f'{pdb_id} assembly {assembly}')
        run(session, f'close #{au.id_string}')

    run(session, f'delete solvent & #{model.id_string}')
    run(session, f'delete (~protein & ~nucleic) & #{model.id_string}')
    extra = spec.get('delete')
    if extra:
        run(session, f'delete ({extra}) & #{model.id_string}')

    # Add hydrogens if the structure is essentially unprotonated. ISOLDE's AMBER
    # templates require explicit H; sub-Angstrom structures (1EJG/2VB1) ship with
    # them, but most don't. This is exactly the recovery ISOLDE itself offers on a
    # bad-template hit (isolde.py::_handle_bad_template -> cmd_addh(hbond=True)).
    # Guarded so already-protonated structures are left untouched (keeps them
    # reproducible and avoids double-adding).
    n_h = int((model.atoms.element_names == 'H').sum())
    n_heavy = model.num_atoms - n_h
    if n_h < 0.5 * max(n_heavy, 1):
        from chimerax.atomic import AtomicStructures
        from chimerax.addh import cmd as addh_cmd
        addh_cmd.cmd_addh(session, AtomicStructures([model]), hbond=True)

    # TODO: alternate-conformer handling. The curated structures are chosen partly
    # for having little/no alt-loc disorder, but a more thorough pass should
    # explicitly pick one conformer (ChimeraX's `altlocs` command) rather than
    # relying on the source structure being clean.
    return model


def _minimize_and_get_virial(session, model):
    '''Build a whole-model (nothing fixed) ISOLDE simulation, minimize once at
    ISOLDE's own default tolerance, and return (atoms, symmetrised_virial).'''
    from chimerax.atomic import Atoms
    # Absolute imports: tools/ is a dev-only directory, not a submodule of the
    # installed chimerax.isolde package (cf. tools/gen_rs_glyphs.py), so relative
    # imports from here don't resolve -- this script must run inside ChimeraX
    # with the isolde bundle already installed and importable.
    from chimerax.isolde.openmm import openmm_interface, sim_param_mgr

    all_atoms = model.atoms
    empty = Atoms()
    sim_construct = openmm_interface.SimConstruct(model, all_atoms, empty, None)

    sim_params = sim_param_mgr.SimParams()
    # 'CPU' is the multithreaded, always-available OpenMM platform -- much faster
    # than 'Reference' for the larger structures here, and platform choice doesn't
    # materially affect a minimised-geometry median baseline (tiny numerical
    # differences well below the stress signal of interest).
    sim_params.platform = 'CPU'

    # session.isolde lazily creates/returns ISOLDE's session singleton, which
    # owns the cached ForcefieldMgr used by every real simulation
    # (openmm_interface.py:329-330) -- reuse it rather than building a second one.
    forcefield_mgr = session.isolde.forcefield_mgr

    sim_handler = openmm_interface.SimHandler(session, sim_params, sim_construct, forcefield_mgr)
    # Deliberately skip initialize_restraint_forces(): with every restraint type
    # off, no restraint Force objects are added at all, and CMAP off keeps the
    # geometry consistent with what LocalVirialCalculator actually reconstructs
    # (it doesn't cover CMAP -- see local_virial.py's module docstring).
    sim_handler.initialize_restraint_forces(
        amber_cmap=False, tugging=False, position_restraints=False,
        distance_restraints=False, adaptive_distance_restraints=False,
        dihedral_restraints=False, adaptive_dihedral_restraints=False)

    # _prepare_sim() is exactly start_sim()'s first step (openmm_interface.py:1549)
    # -- builds the Context (applying the softcore nonbonded conversion, since
    # that's gated on params.use_softcore_nonbonded_potential, not on anything
    # above) -- but skips the rest of start_sim()'s async, frame-triggered
    # continuation, which this one-shot headless pass doesn't want.
    sim_handler._prepare_sim()
    th = sim_handler._thread_handler
    th.minimize()  # ISOLDE's own default tolerance/max_iterations (params-derived)
    th.finalize_thread()  # block until the background C++ minimisation thread is done

    from chimerax.isolde.validation.local_virial import LocalVirialCalculator
    state = sim_handler._context.getState(getPositions=True)
    calc = LocalVirialCalculator(sim_handler._system)
    virial = calc.compute(state.getPositions(asNumpy=True))['virial']
    virial_sym = 0.5 * (virial + np.transpose(virial, (0, 2, 1)))
    return sim_construct.all_atoms, virial_sym


def build(session, structures=None, out_path=None):
    '''
    Build and save a ReferenceStressLibrary cache from a curated structure set.

    Args:
        session: the ChimeraX session (with ISOLDE already available as
            session.isolde).
        structures: iterable of structure specs. Each may be a bare PDB-ID string
            or a dict (see DEFAULT_STRUCTURES for the accepted keys: 'id',
            'assembly', 'delete'). Defaults to DEFAULT_STRUCTURES.
        out_path: where to write the .npz cache. Defaults to
            validation/data/reference_virial.npz alongside this tools/ directory.
    '''
    from chimerax.core.commands import run
    from chimerax.isolde.validation.reference_stress import ReferenceStressLibrary

    if structures is None:
        structures = DEFAULT_STRUCTURES
    if out_path is None:
        import os
        out_path = os.path.join(
            os.path.dirname(__file__), '..', 'src', 'validation', 'data',
            'reference_virial.npz')

    samples = {}
    for spec in structures:
        spec = _normalize_spec(spec)
        pdb_id = spec['id']
        session.logger.info(f'build_reference_stress_library: processing {pdb_id}...')
        model = None
        try:
            model = _prepare_structure(session, spec)
            atoms, virial_sym = _minimize_and_get_virial(session, model)
            n_added = 0
            for i, atom in enumerate(atoms):
                if _is_terminal_residue(atom.residue):
                    continue
                key = f'{atom.residue.name}:{atom.name}'
                samples.setdefault(key, []).append(virial_sym[i])
                n_added += 1
            session.logger.info(
                f'build_reference_stress_library: {pdb_id} contributed {n_added} '
                'atom samples.')
        except Exception as e:
            detail = str(e)
            for attr in ('unmatched', 'ambiguous'):
                residues = getattr(e, attr, None)
                if residues:
                    detail += (f'; {attr}=' + ','.join(
                        f'{r.name}{r.number}/{r.chain_id}' for r in residues[:40]))
            session.logger.warning(
                f'build_reference_stress_library: skipping {pdb_id} ({detail}).')
        finally:
            if model is not None and not model.deleted:
                run(session, f'close #{model.id_string}')

    lib = ReferenceStressLibrary.from_samples(samples)
    lib.save(out_path)
    session.logger.info(
        f'build_reference_stress_library: wrote {len(samples)} (residue, atom) '
        f'baselines to {out_path}.')
    return lib
