
_fetched_templates = set()

def find_incorrect_residues(session, model, heavy_atoms_only = True):
    from chimerax.atomic import Residue, mmcif, TmplResidue
    residues = model.residues
    questionable = []
    for r in residues:
        pt = r.polymer_type
        tmpl = None
        if pt in (Residue.PT_AMINO, Residue.PT_NUCLEIC):
            start=True
            end=True
            if pt == Residue.PT_AMINO:
                fa = r.find_atom('N')
                la = r.find_atom('C')
            else:
                fa = r.find_atom('P')
                la = r.find_atom("O3'")
            if fa is not None:
                for fn in fa.neighbors:
                    if fn.residue != r:
                        start=False
                        break
            if la is not None:
                for ln in la.neighbors:
                    if ln.residue != r:
                        end=False
                        break
            try:
                tmpl = TmplResidue.get_template(r.name, start=start, end=end)
            except ValueError:
                tmpl = None
        if tmpl is None:
            if r.name not in _fetched_templates:
                session.logger.info('Fetching CCD definition for residue {} {}{}'.format(
                    r.name, r.chain_id, r.number
                ))
            try:
                tmpl = mmcif.find_template_residue(session, r.name)
                _fetched_templates.add(r.name)
            except ValueError:
                session.logger.warning('Template {} not found in the Chemical Components Dictionary'.format(r.name))
                continue
        if heavy_atoms_only:
            ra_names = set(r.atoms[r.atoms.element_names != 'H'].names)
            ta_names = set([a.name for a in tmpl.atoms if a.element.name != 'H'])
        else:
            ra_names = set(r.atoms.names)
            ta_names = set([a.name for a in t.atoms])
        ra_residuals = ra_names.difference(ta_names)
        ta_residuals = ta_names.difference(ra_names)
        if len(ta_residuals):
            if end:
                if pt == Residue.PT_AMINO:
                    print('C-terminal residue {} {}{}; ra_residuals: {}; ta_residuals: {}'.format(
                        r.name, r.chain_id, r.number, ra_residuals, ta_residuals
                    ))
                if pt == Residue.PT_AMINO and not len(ra_residuals) and ta_residuals == set(('OXT',)):
                    # Dangling C-terminal peptide. Allow.
                    continue
                elif pt == Residue.PT_NUCLEIC and not heavy_atoms_only and not len(ra_residuals) and ta_residuals == set(("HO5'",)):
                    # Dangling 5' end. Allow
                    continue
            if start and not heavy_atoms_only:
                if pt == Residue.PT_AMINO and ta_residuals == set(('H',)) and ra_residuals == set(('H1','H2','H3')):
                    # Dangling N-terminal peptide. Allow
                    continue
                elif pt == Residue.PT_NUCLEIC and ta_residuals == set(("HO3'",)):
                    # Dangling 3' end. Allow.
                    continue

            questionable.append(r)
    return questionable
