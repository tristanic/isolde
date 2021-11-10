from chimerax.core.commands import AtomSpecArg

class IsoldeStructureArg(AtomSpecArg):

    name = "a single atomic structure"

    @classmethod
    def parse(cls, text, session):
        aspec, text, rest = super().parse(text, session)
        models = aspec.evaluate(session).models
        from chimerax.atomic import AtomicStructure
        # Need to use explicit type comparison here rather than isinstance() since some
        # things we *don't* want (e.g. ISOLDE's rotamer preview) are AtomicStructure
        # subclasses.
        mols = [m for m in models if type(m) == AtomicStructure]
        if len(mols) != 1:
            from chimerax.core.commands import AnnotationError
            raise AnnotationError(f'Must specify exactly one atomic structure, got {len(mols)} for {text}.')
        return mols[0], text, rest
        
