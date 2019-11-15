# @Author: Tristan Croll <tic20>
# @Date:   11-Jun-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 14-Jun-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll



def rota(session, structures=None, report=False):
    '''
    Add a live rotamer validator to each of the given structures, and optionally
    report a summary of current outliers.
    '''
    from chimerax.isolde import session_extensions as sx
    if structures is None:
        from chimerax.atomic import AtomicStructure
        structures = [m for m in session.models.list() if type(m)==AtomicStructure]
    for structure in structures:
        sx.get_rota_annotator(structure)
    if report:
        from chimerax.atomic import Residues, concatenate
        residues = concatenate([m.residues for m in structures])
        mgr = sx.get_rotamer_mgr(session)
        rotamers = mgr.get_rotamers(residues)
        report_str = 'NON-FAVOURED ROTAMERS: \n'
        nf, scores = mgr.non_favored_rotamers(rotamers)
        for r, score in zip(nf, scores):
            report_str += '#{:<6} {}:\t{} {} (P={:.4f})\n'.format(
                r.residue.structure.id_string, r.residue.chain_id, r.residue.name,
                r.residue.number, score
                )
        session.logger.info(report_str)

def unrota(session, structures=None):
    '''
    Delete any rotamer annotators associated with the given models.
    '''
    if structures is None:
        from chimerax.atomic import AtomicStructure
        structures = [m for m in session.models.list() if type(m)==AtomicStructure]
    from chimerax.isolde import session_extensions as sx
    for structure in structures:
        ra = sx.get_rota_annotator(structure, create=False)
        if ra is not None:
            session.models.close([ra])

def rama(session, structures=None, show_favored=True, report=False):
    '''
    Add a live Ramachandran validator to each of the given structures, and
    optionally report a summary of current outliers and cis/twisted peptide
    bonds.
    '''
    from chimerax.isolde import session_extensions as sx
    if structures is None:
        from chimerax.atomic import AtomicStructure
        structures = [m for m in session.models.list() if type(m)==AtomicStructure]
    for structure in structures:
        ra = sx.get_rama_annotator(structure)
        ra.hide_favored = not show_favored
    if report:
        from chimerax.atomic import Residues, concatenate
        residues = concatenate([m.residues for m in structures])
        mgr = sx.get_ramachandran_mgr(session)
        report_str = 'RAMACHANDRAN OUTLIERS: \n'
        outliers = mgr.outliers(residues)
        for outlier in outliers:
            report_str +='#{:<6} {}:\t{} {}\n'.format(
                outlier.structure.id_string, outlier.chain_id, outlier.name, outlier.number)
        report_str += '\nCIS PEPTIDE BONDS: \n'
        cispeps = mgr.cis(residues)
        for cis in cispeps:
            report_str +='#{:<6} {}:\t{} {}\n'.format(
                cis.structure.id_string, cis.chain_id, cis.name, cis.number
            )
        report_str += '\nTWISTED PEPTIDE BONDS: \n'
        twisteds = mgr.twisted(residues)
        for twisted, angle in twisteds:
            report_str += '#{:<6} {}:\t{} {} ({:.1f}°)\n'.format(
                twisted.structure.id_string, twisted.chain_id, twisted.name,
                twisted.number, angle
            )
        session.logger.info(report_str)

def unrama(session, structures=None):
    '''
    Delete any Ramachandran annotators associated with the given models.
    '''
    if structures is None:
        from chimerax.atomic import AtomicStructure
        structures = [m for m in session.models.list() if type(m)==AtomicStructure]
    from chimerax.isolde import session_extensions as sx
    for structure in structures:
        ra = sx.get_rama_annotator(structure, create=False)
        if ra is not None:
            session.models.close([ra])

def register_rota(logger):
    from chimerax.core.commands import (
        register, CmdDesc, BoolArg, create_alias
    )
    from chimerax.atomic import StructuresArg
    desc = CmdDesc(
        optional=[
            ('structures', StructuresArg),
            ],
        keyword=[
            ('report', BoolArg),
        ],
        synopsis='Add rotamer validator markup to models and optionally report current outliers'
    )
    register('rota', desc, rota, logger=logger)
    undesc = CmdDesc(
        optional=[('structures', StructuresArg)],
        synopsis='Close the rotamer annotators for the given models (or all if no models given)'
    )
    register('rota stop', undesc, unrota, logger=logger)
    create_alias('~rota', 'rota stop $*', logger=logger)


def register_rama(logger):
    from chimerax.core.commands import (
        register, CmdDesc, BoolArg, create_alias
    )
    from chimerax.atomic import StructuresArg
    desc = CmdDesc(
        optional=[
            ('structures', StructuresArg),
            ],
        keyword=[
            ('show_favored', BoolArg),
            ('report', BoolArg),
        ],
        synopsis='Add Ramachandran validator markup to models and optionally report current outliers'
    )
    register('rama', desc, rama, logger=logger)
    undesc = CmdDesc(
        optional=[('structures', StructuresArg)],
        synopsis='Close the Ramachandran annotators for the given models (or all if no models given)'
    )
    register('rama stop', undesc, unrama, logger=logger)
    create_alias('~rama', 'rama stop $*', logger=logger)
