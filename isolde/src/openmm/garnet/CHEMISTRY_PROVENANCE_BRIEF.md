# Chemistry provenance & incomplete-model handling for the GARNET force field

**A design brief for future implementation.** Status: design only — nothing here is built yet.
Audience: someone extending ISOLDE's interim GARNET force field (`isolde/src/openmm/garnet/`).
Companion to `FORCEFIELD_ON_MODEL_BRIEFING.md` (the "force field on the model" / Track B briefing);
read this after understanding how the interim GARNET path parameterises the whole model and subsets
it (`params_cache.py`, `system_builder.py`, `primitives.py`).

---

## Context — why this is needed

GARNET parameterises from **topology alone**, in a single forward pass, with no template library. That
is its power, but it removes a safety net the AMBER path relied on: *"this residue matches no template →
refuse to simulate it."* With GARNET, almost **anything** will parameterise and run — including
chemistry that is wrong, mis-protonated, or an artefact of an incomplete model — and it will do so
**silently**.

So the burden shifts to us: we must **confirm up front that each component's chemistry is what we
intend**. Concretely, two requirements:

1. **Provenance.** Every simulated component (residue, ligand, ion, covalent unit) must be connected
   back to a **verified source of truth** — a CCD template, or a graph generated from SMILES / drawn by
   the user — with correct understanding of its inter-residue bonds.
2. **Honest handling of incompleteness.** Real models have missing sidechains and backbone breaks. We
   must support these **without inventing spurious chemistry**. The strategy is to parameterise the
   **expected** graph (what the chemistry *would* be if the model were complete) and then **trim** back
   to what is actually modelled — with explicit rules for which cuts are allowed.

This is a large effort, and it includes a **user-facing graphical framework** for setting, inspecting
and modifying a component's chemistry before a simulation starts. The good news: almost all the
low-level machinery already exists in ISOLDE and the sibling bundles (see [Reuse map](#reuse-map)). The
genuinely new work is the **provenance/gate layer**, the **expected-graph → trim bridge**, the
**cut-rule engine**, and **generalising** the inspector UI.

**Locked design decisions** (with TC):
- The verification is a **hard gate with user override**: a GARNET simulation will not start until every
  component is verified against a source of truth or explicitly overridden.
- The inspector UI **generalises Hamish's** `feature/validation-new-section` widget
  (`ui/validation_tab/new_section.py` + its `ChemSearchButton`), rather than being built from scratch.

---

## 1. Principle

Every simulated component must resolve to a **verified source of truth** and carry a **provenance
record** before GARNET types it. The simulation is **gated**: it will not start until every component
is either verified, or carries an explicit user override.

"Verified" is **not** the same as "complete." An incompletely-modelled component can be verified
*against its expected graph* and simulated on the modelled subset — provided the missing part is an
**allowed cut** (§6).

## 2. Source-of-truth resolution (per component)

Resolve each component to a reference graph, in precedence order. This mirrors the fallback chain
`rdkit_bridge.ccd_records()` already uses:

1. **Registered / overridden** — a per-residue override (`Residue.isolde_chemcomp_id`, or a new
   provenance marker) points at a record in the ChemComp store, including user- or SMILES-registered
   compounds added via `chemcomp.register_record`.
2. **CCD by id** — `chimerax.chemcomp.record()` / `lookup()` (offline SQLite). Gives elements, formal
   charges, bond orders, aromatic flags, stereo config, **leaving-atom flags**, and ideal coordinates.
3. **CCD by graph** — when the name does not resolve or disagrees, graph-match against candidate CCD
   templates (`forcefields.find_possible_templates`, `template_utils.match_ff_templates_to_ccd_templates`).
4. **SMILES / user-drawn** — `rdkit_bridge.records_from_smiles(...)`, or the ChemSearch Ketcher editor →
   `records_from_mol` → `register_record`; bound to the component via `isolde_chemcomp_id`.
5. **Unresolved** → the component is *unknown-source* → hard-blocked until the user supplies a source.

## 3. Per-component verification

For each component, build the **expected graph** and the **modelled graph**, then correspond them:

- **Expected graph:** `TmplResidue.get_template(name, start=, end=)` (with `start`/`end` chosen from the
  residue's inter-residue bonds), or the CCD `record()` connection table; a SMILES/registered source
  supplies its own record.
- **Correspondence:** `template_utils.find_maximal_isomorphous_fragment` — name-match when atom names
  are unique, else a bond-order- and chirality-aware FMCS (RDKit), with `chimerax.isolde.graph`/`mcsplit`
  as a graph fallback.
- **Classification** per component: **verified-complete**, **verified-incomplete (allowed trim)**,
  **disallowed cut**, **mismatch** (modelled atoms / bonds / formal charge / stereo disagree with the
  source), or **unknown source**.

Use `template_utils.find_incorrect_residues` (expected-vs-modelled heavy-atom name-set diff, with its
built-in terminus exceptions) as a fast pre-check; escalate to FMCS for the real decision. Check
chirality and formal charge via `chirality.reference_cip_codes` / `as_modelled_configs` and the CCD
formal charges. This ties directly to the charge pipeline: GARNET's partial charges are **conditioned
on the inferred formal charges** and neutralised per connected component, so a formal-charge disagreement
is a first-class verification failure, not a cosmetic one.

## 4. Inter-residue bonds (expected vs modelled)

Reuse `openmm/amberff/covalent.py`:

- `_is_backbone_continuation` — peptide `C–N` and nucleic `O3'–P` between standard polymers are the
  canonical intra-polymer links: **boundaries**, not errors.
- `_inter_residue_bonds` — graph walk of every bond leaving a residue.
- `detect_covalent_unit` / `CovalentUnit` — grow a unit across **non-standard** seams (disulfides,
  glycosidic, isopeptide, ligand attachments, modified backbones); a covalent unit is verified and
  parameterised **as a unit**, against its combined expected graph.

The CCD **leaving-atom flags** (`pdbx_leaving_atom_flag`, surfaced by `chemcomp.record()` and
`chemcomp.chiral_centres().linkage`) predict *where* inter-residue bonds belong, so the framework can
confirm a modelled inter-residue bond is expected.

Inter-residue anomalies are handled gracefully, **not** by blocking:
- A **missing** expected link is just a **chain break** — treat the residue as a terminus/continuation
  (§5), exactly like the backbone-continuation case.
- An **unexpected** extra link triggers **covalent-unit formation** and verification of that unit;
  it blocks only if the resulting unit is itself unresolvable.

## 5. Expected-graph → trim parameterisation (the core new bridge)

Today GARNET's `subset_for` / `SubsetParameters.missing_atoms` are pure **topology-membership** tests
(is this atom in `all_atoms`?). They cannot distinguish "legitimately outside the simulated region" from
"missing because it was never modelled" — the GARNET cache has **no link** to any expected graph. That
bridge is the central new mechanism:

- **Parameterise the expected graph** of each component (complete residue; a broken backbone extended by
  a chemically-sensible continuation/cap; a covalent unit assembled whole), giving GARNET a
  chemically-complete context. Its degree-based typing and per-component charge neutralisation both need
  this — the same rationale as the current whole-model pass.
- **GARNET is only ever handed this complete graph, never the modelled truncation.** It therefore
  cannot mis-type a truncated fragment; correctness comes from typing the *intended* chemistry and then
  trimming — not from guarding what GARNET sees. (Concretely: a sidechain modelled only to Cβ is typed
  as the full residue, so its Cβ receives the parameters it would have in the intact sidechain.)
- **Trim** back to the actually-modelled atoms, exactly as `subset_for` already trims to `all_atoms`:
  terms with an absent atom are simply dropped.
- The **provenance layer keys expected-graph atoms to modelled ChimeraX atoms** (the missing link), so a
  dropped term is known to be "expected but not modelled" rather than "outside the simulated region."
- The **frozen buffer edge is unchanged** (whole-residue shell + `intra_bonds`-only topology +
  `ignoreExternalBonds`) — that is a separate, legitimate freezing mechanism.

This expected-graph → trim step **obviates ISOLDE's special truncated-sidechain templates**: typing the
full residue and dropping the terms to missing atoms replaces them.

## 6. Cut rules (which trims are allowed)

Because parameterisation is always done on the **expected (complete) graph** (§5) — GARNET never sees
the truncated topology — the cut rules are **not** about what GARNET would make of the modelled atoms.
They turn on two questions: **(i)** is the **expected graph unambiguously determinable** for the
component, and **(ii)** does **trimming to the modelled atoms leave a physically coherent fragment**
(kept atoms remain adequately bonded/restrained; no charged group left half-modelled)?

**Allowed** (parameterise-complete, then trim):
- Chain-terminus differences (OXT, terminal-H variants).
- A **sidechain truncated toward the backbone** (e.g. to Cβ) — verified against the full-residue
  template and typed complete, then trimmed. Terms to the absent atoms are simply **dropped** (the
  existing `subset_for` trim), with **no cap, no filler H, and no forced freezing** — the truncation
  point stays **mobile** unless it independently falls in the fixed shell. This is a common must-support
  case (ChimeraX `addh` deliberately does *not* cap these).
- Generally, any **terminal / pendant** truncation that leaves a connected subgraph of the source.

**Hard-no** (block; require a fix or an explicit override):
- The expected graph is **ambiguous** — an unknown component, or missing atoms that make the intended
  protonation / tautomer / formal charge indeterminate (which also poisons the formal-charge-conditioned
  charge prediction).
- An **interior** hole rather than a terminal truncation — above all a **partially-modelled ring** (the
  surviving ring atoms lose the closure that fixed their geometry and are left under-restrained), or a
  gap inside a conjugated / aromatic system.

*(Inter-residue bonds are deliberately not in this list — see §4: a missing expected link is a chain
break, an unexpected link forms a covalent unit.)*

**Rule of thumb.** A cut is allowed iff **(a)** the complete expected graph is unambiguous, and **(b)**
the kept atoms form a **connected, terminal** subgraph of it (no interior or ring atoms missing) **and
no strongly-charged group is left partially modelled** — a charged group (ammonium, carboxylate,
guanidinium, phosphate, …) must be kept whole or dropped whole. So LYS/ASP/GLU truncated to Cα/Cβ is
allowed (the *entire* charged group is gone), while a half-modelled guanidinium is not.

> This is deliberately **not** the stronger rule "formal charge on the kept atoms is unchanged" — that
> would wrongly forbid the everyday LYS/ASP/GLU truncations.

> **⚠ Open — the phosphodiester carve-out.** Standard nucleic-acid CCD templates place the (charged)
> phosphate *across* the inter-residue boundary, so a naive "charged group all-or-nothing" rule collides
> with normal backbone-continuation handling. The rule needs a boundary-aware exception — a
> backbone-continuation-straddling charged group must be treated differently from a sidechain one. Not
> yet resolved; the trickiest cut-rule corner.

## 7. Gate integration (hard gate + override)

Run a **preflight verification** at simulation construction, mirroring the AMBER path's
`UnparameterisedResiduesError` → widget flow and the read-only `isolde preflight parameters` report. On
any unverified / disallowed-cut component, **refuse to start** and open the inspector.

**Override** is explicit and persisted per residue via session-serialised custom attributes
(`isolde_chemcomp_id`, `isolde_template_name`, plus a new `isolde_ff_provenance` marker) — the same
hook pattern ISOLDE already uses for template overrides.

## 8. UI — generalise Hamish's `new_section.py`

Build on `feature/validation-new-section` (coordinate its merge). That widget already provides:
per-component rows; candidate boxes colour-scored by RDKit FMCS overlap; hover → live 3D preview
(thin-stick); click → rebuild/commit; a `ParameteriseRow` for covalent-unit / metal-site / novel-ligand
builds; and a `ChemSearchButton` (✎) that opens the **ChemSearch Ketcher editor**.

Generalise it from "AMBER MD-template match" to "chemistry provenance":
- A per-component **status** column (verified / incomplete-allowed / disallowed / mismatch / unknown)
  using `QLedLabel`.
- Show the **expected graph with modelled atoms highlighted and missing atoms greyed** — reuse ChemSearch
  `render.py` greying and `chem.svg_data_uri(grey_atoms=…)`, built from `chem.mol_from_ccd` /
  `mol_from_smiles` / `mol_from_geometry`. (ChemSearch already does exactly this "whole molecule, part
  greyed" depiction for contained-substructure hits.)
- The ✎ editor is the **draw / override surface**: SMILES or hand-drawn structure → `records_from_mol`
  → `register_record` → bound to the component via `isolde_chemcomp_id`.
- Every action **logs an `isolde …` command** (the GUI-logs-its-command convention).
- Per-component 3D markup as a `ChiralAnnotator`-style child `Model`.

It lives as a **Validation-tab section**; a hard-block at sim start surfaces it.

## 9. Track B alignment

Provenance and force-field parameters both key to **ChimeraX atom / residue identity**, so this layer is
the natural companion to the future C++ `FFAtom` / `FFBond` / `FFAngle` store described in
`FORCEFIELD_ON_MODEL_BRIEFING.md`. Design the provenance record so it can **migrate onto** a C++ `FF*` /
per-residue object later. For now it lives as Python (beside the `GarnetParameters` cache) plus
session-serialised custom attributes. The expected-graph → trim step is **upstream** of parameter
storage in both the interim and the Track B designs.

## 10. Staged plan

- **P1 — provenance backend (headless-testable).** Source resolution + expected/modelled correspondence
  + classification + cut-rule engine, returning a per-component report. Reuses `ccd_records`,
  `find_maximal_isomorphous_fragment`, `covalent.py`, `chirality`. Verifiable on loaded models without a
  simulation, like `isolde preflight parameters`.
- **P2 — expected-graph → trim parameterisation.** Build the expected graph per component/unit, key it
  to modelled atoms, feed GARNET, trim via the `subset_for` pattern; bridge cache ↔ expected graph.
- **P3 — hard gate.** Preflight at sim start; block + override; persist overrides.
- **P4 — UI.** Generalise `new_section.py`: status list + greyed expected-graph depiction + ChemSearch
  ✎ editor + 3D markup.
- **P5 — Track B migration hook.** Shape the provenance record so it can move onto the C++ FF store
  without a second rewrite.

## 11. Reuse map

- **CCD source of truth:** `chimerax-chemcomp/src/collections.py` — `record`, `lookup`, `chiral_centres`;
  store in `store.py` (offline SQLite; msgpack blobs).
- **Verify vs template:** `isolde/src/atomic/rdkit_bridge.py` (`ccd_records`, `residue_to_rdkit`,
  `template_to_rdkit`, `records_from_smiles`, `records_from_mol`); `atomic/template_utils.py`
  (`find_maximal_isomorphous_fragment`, `find_incorrect_residues`,
  `fix/trim/add_missing_md_template_atoms`, `match_ff_templates_to_ccd_templates`);
  `atomic/chirality.py`; `graph/` + `mcsplit`.
- **Expected graph:** ChimeraX `TmplResidue` / `TmplAtom` (`atomic_cpp/cytmpl.pyx`);
  `mmcif.find_template_residue`.
- **Inter-residue bonds:** `isolde/src/openmm/amberff/covalent.py`.
- **Trim primitive:** `isolde/src/openmm/garnet/params_cache.py` — `subset_for`,
  `SubsetParameters.missing_atoms`.
- **Gate precedent:** `openmm/openmm_interface.py` — `UnparameterisedResiduesError`,
  `_build_system_with_template_recovery`; `validation/cmd.py` — `isolde preflight parameters`.
- **UI:** Hamish `feature/validation-new-section` — `ui/validation_tab/new_section.py` (+ `ChemSearchButton`);
  ChemSearch `chem.py` / `render.py` / `tool.py` (Ketcher editor, atom greying); ISOLDE `ui/ui_base.py`
  (`UI_Panel_Base`, `QLedLabel`, `ExpertModeSelector`), `ui/validation_tab/chiral.py` (per-row
  action-button + click-to-focus pattern), `validation/chiral_annotation.py` (child-`Model` 3D markup).

## 12. Open questions & risks

- **Backbone-break continuations.** What is the "expected" extension at a chain break — a capped standard
  residue, a hydrogen cap, or a user choice? Needs a sensible default plus an override.
- **Charged-group cut rule vs the phosphodiester convention.** "A strongly-charged group is cut whole"
  handles LYS/ASP/GLU, but nucleic-acid CCD templates split the charged phosphate across the residue
  boundary (`O3'–P`), so the rule needs a boundary-aware exception. Design TBD; the trickiest corner.
- **SMILES / registered-source stereo & protonation** must be pinned down (ChemSearch already carries
  dimorphite-dl pH-7 handling to draw on).
- **Coordinating the merge** of Hamish's branch vs generalising on top of it.
- **Performance.** FMCS is NP-hard; Hamish's widget already caps per-comparison cost. Memoise CCD
  lookups per id.
- **Track B migration** of the provenance record onto the C++ store without a second rewrite.

## 13. Verification of the future implementation

- **P1** is headless-testable on loaded models (no simulation), like `isolde preflight parameters`:
  assert the per-component classification and cut-rule verdicts on curated cases (complete residue;
  Cβ-truncated LYS/ASP/GLU → *allowed*; partially-modelled ring → *hard-no*; unknown ligand → *unknown
  source*; a disulfide / glycosidic link → *covalent unit*).
- **P2** cross-checks that a component parameterised on its expected graph and trimmed matches, on the
  kept atoms, the parameters it gets when the component is modelled complete.
- **P3–P4** are GUI-only (per ISOLDE's no-headless-sim rule): the gate blocks and the inspector opens on
  an unverified component; override lets the sim proceed and persists across save/reload; the greyed
  expected-graph depiction and the ✎ editor round-trip a drawn/registered override.
