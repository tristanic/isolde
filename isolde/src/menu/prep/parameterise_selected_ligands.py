# @Author: Tristan Croll <tic20>
# @Date:   04-Aug-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 04-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

tooltip = ('Parameterise the selected residues for MD. Hydrogens must be present and correct.')

def run_script(session):
    from chimerax.core.commands import run
    run(session, 'isolde parameterise sel')
