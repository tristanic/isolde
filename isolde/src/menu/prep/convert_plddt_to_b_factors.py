tooltip = ('<span>Convert values in the B-factor column of the working model from pLDDT to '
    'estimated B-factors. This should <strong>only</strong> be used if the model is '
    'an AlphaFold prediction.</span>')

display_name = 'Convert pLDDTs to B Factors'

def run_script(session):
    from chimerax.core.errors import UserError
    from chimerax.isolde.reference_model.alphafold import convert_plddt_to_b_factors
    isolde = getattr(session, 'isolde', None)
    if isolde is None:
        raise UserError('ISOLDE must be running first!')
    m = isolde.selected_model
    if m is None:
        raise UserError('ISOLDE does not have a model currently selected!')
    convert_plddt_to_b_factors(m)
    session.logger.info(f'ISOLDE: Converted pLDDT values to B-factors for model #{m.id_string}')
