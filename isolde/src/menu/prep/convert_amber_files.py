# @Author: Tristan Croll <tic20>
# @Date:   15-Jul-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 03-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

tooltip = ('Convert all pairs of .frcmod/.mol2 AMBER definition files in the current directory to the OpenMM ffXML format')

def run_script(session):
    print('running convert_amber_files')
    from chimerax.isolde.parmed import install_parmed_if_necessary
    install_parmed_if_necessary(session)
    from chimerax.isolde.openmm.amberff.amber_convert import amber_to_xml_individual
    amber_to_xml_individual('.','.', search_subdirs=False)
