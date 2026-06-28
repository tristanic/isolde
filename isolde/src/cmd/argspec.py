from chimerax.core.commands import AtomSpecArg, AnnotationError


def isolde_target_structures(models):
    '''
    Filter *models* down to the genuine atomic structures ISOLDE should act on,
    excluding the ``AtomicStructure``-derived helper models that live *under* a
    working structure and must never be validated/annotated/simulated:

      * ISOLDE's rotamer preview -- a ``_RotamerPreview(AtomicStructure)``
        **subclass**; and
      * ChimeraX's altloc-display drawings -- plain ``AtomicStructure`` instances
        nested in an "alternate locations" group beneath the model (added by
        recent ChimeraX versions; ``altlocs show``).

    The first is excluded by the exact-type test (it is a subclass, not
    ``AtomicStructure`` itself). The second is excluded the same way ChimeraX's
    own ``match_maker._remove_child_models`` does it: a candidate that appears
    inside *another* candidate's model tree (``all_models()``) is a sub-structure,
    not a working model. This is safe for the real working model -- it can only be
    dropped if it sits under another candidate structure, which never happens.
    '''
    from chimerax.atomic import AtomicStructure
    # Exact type, not isinstance(): subordinate helpers like _RotamerPreview are
    # AtomicStructure *subclasses* we want to drop.
    candidates = [m for m in models if type(m) is AtomicStructure]
    cset = set(candidates)
    nested = set()
    for s in candidates:
        for child in s.all_models():
            if child is not s and child in cset:
                nested.add(child)
    return [s for s in candidates if s not in nested]


class IsoldeStructureArg(AtomSpecArg):

    name = "a single atomic structure"

    @classmethod
    def parse(cls, text, session):
        aspec, text, rest = super().parse(text, session)
        mols = isolde_target_structures(aspec.evaluate(session).models)
        if len(mols) != 1:
            raise AnnotationError(
                f'Must specify exactly one atomic structure, got {len(mols)} for {text}.')
        return mols[0], text, rest


class IsoldeStructuresArg(AtomSpecArg):

    name = "one or more atomic structures"

    @classmethod
    def parse(cls, text, session):
        aspec, text, rest = super().parse(text, session)
        mols = isolde_target_structures(aspec.evaluate(session).models)
        if not mols:
            raise AnnotationError(f'No atomic structures matched {text}.')
        return mols, text, rest
