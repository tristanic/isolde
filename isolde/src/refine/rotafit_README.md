# `rotafit` — settle-and-rank fitting into density, and the ISOLDE-engine playbook behind it

This document covers **`rotafit.py`**: what it does, the algorithm, and — most importantly
for future work — the **hard-won lessons about driving ISOLDE's live OpenMM simulation from
a command**. `rotafit` is the first of a family of "place a candidate, let the force field
+ density decide" tools; the next members (**building residues along a chain**, **docking
ligands with or without covalent bonds**) will reuse almost all of the machinery and every
one of the engine lessons here. Read the "Engine playbook" section before writing any of
them — most of it was learned the expensive way.

---

## 1. What `rotafit` does and why

**Command:** `isolde rotafit [residues] [options]` (single residue by default; see
`allowMultiple`).

**Problem it solves.** ISOLDE's backrub rotamer fitter rotates the sidechain about a *fixed*
CA–CB bond and lets the human eyeball which preview best matches the density. That fails on
a recurring class of cases (e.g. a valine or isoleucine fitted *backwards* into a low-res
map) because:
- the rigid preview can put even the *correct* rotamer well out of the map (the CA–CB bond
  isn't where it should be until the whole thing relaxes), so eyeballing pre-settle previews
  is misleading; and
- long rotamer lists (Arg, Lys) are sorted by *prevalence*, not similarity to the density.

**Core idea.** Don't judge a rotamer by how it looks *before* relaxation — judge it by its
energy *after* the full ISOLDE force field (bonded + soft-core nonbonded + MDFF map) has
settled it. Once a rotamer is placed in roughly the right basin, the force field reliably
falls into the correct conformation, so the **settled energy is a far more honest arbiter
than the pre-settle appearance**. `rotafit` automates the tedious middle step of the manual
workflow (start sim → dial through previews → commit the best).

---

## 2. The algorithm (pipeline)

Per target residue, inside the **running, paused** simulation:

1. **Enumerate candidates.** Every library rotamer (via the rotamer manager, the same
   conformations the preview buttons cycle) **plus the residue's *current* conformation**.
   The current conformation is a first-class candidate (see "do no harm").
2. **Cull non-starters.** Drop any rotamer whose moved heavy atoms severely overlap a fixed
   environment heavy atom (`SEVERE_OVERLAP`, cheap `find_close_points` pre-filter). The
   current conformation is *exempt* from culling.
3. **Soft search.** Settle each survivor at a softened soft-core λ (`settleLambda`, default
   0.6) so a rotamer seeded into an overlap slides apart instead of exploding. Rank by
   scoring energy (see §4).
4. **Polish + re-rank.** Take the top `polishTop` survivors (plus always the current), and
   **ramp** the soft-core λ from `settleLambda` up to full stiffness over `rampIncrements`
   stages (adiabatic — a sudden jump jolts atoms in a soft overlap that becomes a hard
   wall). Re-rank by the polished energy. The soft and stiff rankings can disagree; the
   stiff one is the honest one.
5. **"First, do no harm."** Keep the current conformation unless a rotamer beats it by a
   real margin (`acceptMargin`, scaled by the MDFF coupling — see §4). This stops settling
   noise from flipping an already-correct fit on repeat calls.
6. **Commit (graft).** Commit **only the target residue** onto the *original* environment
   (see §3.7). Leave the simulation running so the live MDFF trajectory refines further and
   the user can `isolde sim revert` to the checkpoint taken at the start.

Settling uses a **hybrid** relax (`minimize=True`, now the default): step the integrator at
0 K first (dynamics has momentum to seat the pose in the local density well), then
`minimizeEnergy()` to converge deterministically (kills the run-to-run jitter that makes
the ranking energy noisy). See §3.6.

---

## 3. Engine playbook — driving ISOLDE's live sim from a command

This is the valuable part. ISOLDE runs OpenMM on a **worker thread** behind a
`SimHandler`/`OpenmmThreadHandler`; a command runs on the **main thread**. Every rule below
was a bug we hit.

### 3.1 Pausing is asynchronous — defer to the `'sim paused'` trigger
`sim_handler.pause = True` only *sets a flag*; the worker keeps running until it next
returns. Driving the simulation immediately after races the worker. **Do the work in a
one-shot handler on the `'sim paused'` trigger** (or immediately if already paused).

Do **not** use the `'coord update'` trigger: it fires *inside*
`_update_coordinates_and_repeat` (openmm_interface.py) *before* that function's pause/repeat
branch. Resuming from there makes the pause-setter schedule one step-loop and the still-live
frame schedule another — two overlapping loops leak per call, each firing `'sim terminated'`
on shutdown (the `remove_change_cb` traceback cascade). `'sim paused'` fires from a
`'frame drawn'` handler *after* that function has returned down the pause branch, so the sim
is genuinely quiescent and resuming schedules exactly one loop.

### 3.2 Drive the integrator directly on the main thread
In the paused window, step/minimise the OpenMM context **directly** — `sim_handler._main_integrator.step(n)`
and `sim_handler._simulation.minimizeEnergy()`. Do **not** use `thread_handler.step/minimize`
(those are threaded/async — more careful code, and unnecessary here). `_main_integrator` is a
plain `LangevinIntegrator` (no `minimizeEnergy` of its own); minimise via `_simulation`.
`_context` **is** `_simulation.context` (one context; no dual-context subtlety).

Stepping `_main_integrator` directly also **bypasses ISOLDE's fast-atom surveillance
integrator** in the CompoundIntegrator — which is what we want: it would otherwise veto/
interrupt a deliberately clashy trial before we score and discard it.

### 3.3 Pushing coordinates while paused
`sim_handler.push_coords_to_sim(coords)` only *queues* the update when paused. Follow it
immediately with the private `sim_handler._push_coords_to_sim()` to flush it to the
thread-handler's canonical buffer so it survives resume. (Helper: `_push_now`.) Coordinates
are in **Ångström**; read them back with `getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.angstrom)`.

### 3.4 Temperature and soft-core λ
- Temperature: `sim_handler.temperature = value` (routes to `_main_integrator`; the
  CompoundIntegrator has no get/setTemperature). Read the *configured* value to restore
  from `isolde.sim_params.temperature`.
- Softening: `rotafit` uses **per-group soft-core coupling** (`SimHandler.assign_nb_group` /
  `set_nb_coupling`, via the `_Softener` helper), NOT the global `softcore_lambda`. The target
  residue goes in its own nonbonded group and only its coupling to the rest of the model is
  softened (~0.6) so a clashy seed relaxes, while the environment stays rigid at full strength
  (no position pinning). **Always restore** temperature and the coupling in a `finally`. (The
  global `sim_handler.softcore_lambda` still acts as the master ceiling; see §3.8.)

### 3.5 Scoring on force groups — exclude restraints
`getState(getEnergy=True).getPotentialEnergy()` returns the **total**, which includes
ISOLDE's `RESTRAINT_FORCE_GROUP` (=6) — and *every* restraint lands there (omega / phi-psi /
distance / chiral / …), pose-independent **noise** for this decision. Rank on
**`CORE_FORCE_GROUPS | {MAP_FORCE_GROUP}` = groups {0,1,2,3,4,5}** (bonded + nonbonded +
MDFF), never group 6. `getState(getEnergy=True, groups=set_of_ints)` accepts a set. Force
groups are defined at the top of openmm_interface.py.

### 3.6 Settling: `minimize` under-seats; use the hybrid
A pure `minimizeEnergy()` from a freshly chi-rotated pose is a **local 0 K quench with no
momentum**: it converges to the *nearest* minimum, which under the gentle cryo-EM map
(`global_k ~ 0.34`) is often **off-density** — it can't ride the small barrier into the
density well. Dynamics *can* (momentum). So the hybrid: **step first (dynamics seats the
pose in the density well), then minimise (converges deterministically to *that* seated
minimum).** This gives clean, noise-free ranking energies *and* density-seated geometry.
Pure-dynamics (`minimize=False`) seats well but its stopping point is jittery (bad for
ranking); pure-minimise is clean but off-density (bad geometry). Chain them.

### 3.7 Committing: graft the target onto the original environment
The settled/polished poses are **whole-construct snapshots** whose *environment* was
perturbed (settled/minimised) during evaluation. Committing that wholesale moves
the surroundings too — invisible for a short 0 K step, but a minimise shifts them enough to
look like a bad commit. **Commit only the target residue, grafted onto the original
coordinates (`base`):** `commit = base.copy(); commit[ridx] = best[ridx]`. Keeping current =>
commit `base` unchanged (literally no change — true do-no-harm). The live sim resolves any
interface afterward.

**Corollary — don't infer per-residue fit from the whole-construct map term.** The map
energy is dominated by the hundreds of in-density environment atoms; a −10000-ish map term
says nothing about whether the single target residue is in density.

### 3.8 Softening the target: per-group soft-core coupling (no pinning)
`rotafit` puts the target residue's real atoms in nonbonded group 1 and its symmetry copies
in group 2 (`SimHandler.assign_nb_group(atoms, 1, copy_group_id=2)`), leaving everything else
in group 0. It then softens every coupling that **touches** the target — vs environment
`(1,0)`/`(2,0)` and its crystal self-contact `(1,2)`/`(2,2)` — while the target's **internal**
`(1,1)` and the environment's `(0,0)` stay **full**. So the sidechain relaxes into density
(soft against its surroundings, *including its own symmetry image* across an interface) while
the environment holds its shape by its own force field — **no position restraints**. Couplings
update live on the paused context (no reinit); the polish ramps them back to full. The
`_Softener` helper wraps this (`assign_target` once, `set(lam)` ramped repeatedly); the
one-call equivalent is `SimHandler.soften_nb_selection`. Needs `nb_groups_max >= 3` (default 4).

> Why per-group and not the global λ: softening the global `softcore_lambda` softens
> *everything*, which used to force `rotafit` to pin the surrounding model with position
> restraints so it wouldn't deform under the soft λ. Per-group coupling never softens the
> environment in the first place, so the pinning machinery (and its risk of *preventing*
> legitimate environmental accommodation, §5) is gone.

### 3.9 Coupling-scaled thresholds
`acceptMargin` is a **multiple of the MDFF coupling constant** (summed `global_k` over
`isolde.sim_manager.mdff_mgrs`), not an absolute energy. `global_k` is sigma-normalised and
resolution-calibrated (small for high-sigma cryo-EM, larger for x-ray), so a coupling-scaled
threshold tracks the map's energy scale automatically across map types. Any energy threshold
compared against MDFF-influenced scores should be scaled this way, not hard-coded in kJ/mol.

### 3.10 Numpy-in-tuples gotcha
Result tuples carry numpy coordinate arrays: `(energy, name, coords, is_current)`. Never use
`x in list` or `==` on them (ambiguous truth value). Track membership with an explicit bool
flag and `any(t[3] for t in ...)` / `next((t for t in ... if t[3]), None)`.

---

## 4. Scoring philosophy

- Rank on **core FF + MDFF** (§3.5); the map influence is ISOLDE's own coupling — there is
  **no separate map weight to invent**. This was a temptation worth resisting: the sim-based
  scorer already inherits the tuned, sigma-normalised coupling.
- The score is a **whole-construct** energy, but the *signal* is one residue. The dominant
  risk is that environment jitter (mostly in the O(N²) nonbonded term) swamps the tiny
  per-residue signal. Mitigations in place: the environment stays rigid at full coupling
  (§3.8), and the deterministic minimise (§3.6) removes stochastic jitter.
- Soft (search) and stiff (polish) rankings can disagree; polish + re-rank at full stiffness
  is the honest call, and the ramp avoids the jolt.

---

## 5. Known limitations

- **Environmental accommodation.** `rotafit` evaluates each rotamer against the environment
  *as it currently sits*. When the correct fix requires the *neighbours* to collectively
  rearrange (cross barriers to a new packing) — e.g. a backwards-fit VAL whose surroundings
  have adapted to the wrong pose — a quick settle can't discover it: the correct pose reads
  as strained against the un-relaxed environment. Giving the sim more room
  (ISOLDE's standard, non-contracted sim-start selection — now the default) helps; a fully
  contracted/over-constrained region can leave no room at all. But some cases genuinely need the
  human to place + equilibrate, with `rotafit` then confirming.
- **Density can't always discriminate.** Near-symmetric residues (VAL, LEU, THR) at low
  resolution put almost identical density for "correct" vs "backwards"; the per-residue map
  signal is near-silent, so the environment (and its incumbency bias) dominates.
- It is a **basin-finder + advisory committer**, human-triggered per residue — not an
  apply-to-all refiner. The `allowMultiple` guard enforces this by default.

---

## 6. Guidance for the next efforts (residue building, ligand docking)

The pipeline generalises to **"place a (possibly tethered) fragment, let FF + density
decide."** What transfers directly:

- **The whole engine playbook (§3).** Threading, direct-integrator driving, coord push,
  force-group scoring, the seat-then-converge hybrid, graft-on-commit, coupling-scaled
  thresholds, and **per-group soft-core coupling** (soften a fragment against its surroundings
  — including its symmetry copies — while keeping the environment rigid, via
  `SimHandler.soften_nb_selection` / the `_Softener` pattern; §3.8). These are engine facts,
  not rotafit-specific — reuse them verbatim.
- **Enumerate → cull → soft-search → polish → do-no-harm → graft-commit** is a generic loop
  over *candidate placements*. For rotamers the candidates are chi-rotations; for a **new
  residue** they're backbone/rotamer placements consistent with the chain geometry; for a
  **ligand** they're docked poses (and conformers). Only the candidate *generator* changes.
- **Do-no-harm** matters even more for building/docking: always include "leave it / don't
  build" as a candidate so a bad placement isn't forced.
- **Graft-on-commit** is essential for docking: commit only the ligand/new-residue atoms onto
  the untouched environment; never commit the evaluation-perturbed neighbourhood.
- **Covalent ligands** are the "tethered fragment" case — like a sidechain, it's bonded to the
  model, so the flanking/junction geometry must be owned by the force field (build the bond
  into the sim topology), not approximated by restraints.

What needs new design:
- **Candidate generation** (poses/conformers) — likely an RDKit-layer job (bond order, ring
  membership, conformer enumeration), dovetailing with the `rdkit_bridge`/`chemcomp` work.
- **Scoring for a fragment that isn't yet in the topology** — a new residue/ligand must be
  added to the simulation construct before it can be settled/scored; decide whether to
  rebuild the context (expensive) or pre-include a placeholder.
- **The environmental-accommodation limitation (§5) will be worse** for bulky ligands — plan
  for "place, then let the pocket relax" (real dynamics), not just a quick settle.

---

## 7. File / symbol map

| Symbol | Role |
|---|---|
| `rotafit()` | command entry: validate, guard (`allowMultiple`), start sim if needed, schedule work on `'sim paused'` |
| `_run_rotafit()` | the worker (runs only when the sim thread is idle): checkpoint, set coupling/temp, per-residue group→settle→commit, restore + resume |
| `_settle_and_rank()` | enumerate + cull + soft-search a residue's candidates; returns ranked results + `base` |
| `_commit_best()` | polish top-N (ramp), re-rank, do-no-harm, **graft-commit target onto `base`** |
| `_settle_pose()` / `_ramped_polish()` | one settle / a λ-ramped polish of a single candidate |
| `_relax()` | the hybrid: 0 K step (seat) then optional `minimizeEnergy()` (converge) |
| `_score_energy()` / `_score_groups()` / `_energy_breakdown()` | force-group–filtered scoring (core FF + MDFF; excl. restraints) + debug breakdown |
| `_Softener` | the softening knob: target real atoms → group 1, copies → group 2; softens every coupling touching the target except its internal (1,1); no pinning |
| `_sim_coupling_constant()` | summed MDFF `global_k` for coupling-scaled `acceptMargin` |
| `_push_now()` / `_sim_coords()` | immediate coord push while paused / read sim coords (Å) |
| `_summary_line()` | the single non-debug log line |

Force-group constants (`CORE_FORCE_GROUPS`, `MAP_FORCE_GROUP`, `RESTRAINT_FORCE_GROUP`) and
the `SimHandler`/`OpenmmThreadHandler` live in `../openmm/openmm_interface.py`.
