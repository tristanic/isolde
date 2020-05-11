# @Author: Tristan Croll <tic20>
# @Date:   06-May-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 06-May-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

def add_isolde_citation(model):
    from chimerax.atomic import mmcif
    mmcif.add_citation(model, 'isolde2018', {
        'title':    'ISOLDE: a physically realistic environment for model building into low-resolution electron density maps',
        'pdbx_database_id_DOI':     '10.1107/S2059798318002425',
        'journal_abbrev':   'Acta Cryst. D',
        'journal_volume':   '74',
        'page_first':   '519',
        'page_last':    '530',
        'journal_issue':    '6',
        'year':     '2018',
        'pdbx_database_id_PubMed':  '29872003'
    })

    from chimerax.isolde import __version__
    import os
    mmcif.add_software(model, 'ISOLDE', {
        'name': 'ISOLDE',
        'version':  __version__,
        'location': 'https://isolde.cimr.cam.ac.uk',
        'classification':   'model building',
        'os':   os.sys.platform.title(),
        'type': 'program',
        'citation_id':  'isolde2018',

    })
