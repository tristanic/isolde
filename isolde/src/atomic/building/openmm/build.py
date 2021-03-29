# @Author: Tristan Croll <tic20>
# @Date:   26-Jun-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 26-Jun-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

def fix_from_openmm_template(residue, md_template, coord_template=None):
    '''
    Add/remove atoms from a residue to make it match the given OpenMM template.
    *NOTE*: MD templates do not provide coordinates nor know anything about
    chirality, so filling in atoms from the MD template alone may yield
    undesirable results. Where possible, provide a ChimeraX TmplResidue object
    as the optional coord_template argument.
    '''
    if coord_template is not None:
        # Fill in all atoms we can from the coord_template first, then
        # make sure it's consistent with the md template.
        from ..template_utils import fix_residue_from_template
        fix_residue_from_template(residue, coord_template)
        return fix_from_openmm_template(residue, md_template)
