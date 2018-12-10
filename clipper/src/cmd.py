
def open_mtz(session, path, structure_model = None,
        over_sampling=1.5, always_raise_errors=False):
    if structure_model is None:
        if always_raise_errors:
            raise TypeError('Reflection data must be associated with an atomic '
                'structure, provided via the structure_model argument.')
        else:
            session.logger.warning('structureModel argument is required when '
                'opening a reflection data file!')
    from .symmetry import get_map_mgr
    mmgr = get_map_mgr(structure_model, create=True)
    try:
        xmapset = mmgr.add_xmapset_from_mtz(path, oversampling_rate=over_sampling)
        log_str = 'Opened crystallographic dataset from {}\n'.format(path)
        if xmapset.experimental_data:
            log_str += 'Found experimental reflection data: \n'
            log_str += '\n'.join(['\t{}'.format(n) for n in xmapset.experimental_data.keys()])
            log_str += '\n'
            log_str += 'Rwork: {:.4f}; Rfree: {:.4f}'.format(
                xmapset.rwork, xmapset.rfree
            )
        log_str += 'Generated maps: \n{}\n'.format(
            '\n'.join(['\t{}'.format(m.name) for m in xmapset]))
        log_str += 'Any unwanted maps may be safely closed via the Model panel.'
        return [mmgr.crystal_mgr], log_str

    except RuntimeError as e:
        if always_raise_errors:
            raise e
        else:
            session.logger.warning(str(e))
            return None, None

def spotlight(session, models=None, enable=True):
    from chimerax.clipper.symmetry import get_symmetry_handler
    if models is None:
        from chimerax.atomic import AtomicStructure
        models = session.models.list(type=AtomicStructure)
    for m in models:
        sh = get_symmetry_handler(m, create=True)
        session.logger.info('Setting spotlight mode for model {} to {}'.format(
            m.id_string, enable
        ))
        sh.spotlight_mode=enable

def associate_volumes(session, volumes, to_model=None):
    if to_model is None:
        from chimerax.core.errors import UserError
        raise UserError('The toModel argument must be provided!')
    from chimerax.clipper.symmetry import get_map_mgr
    mgr = get_map_mgr(to_model, create=True)
    for v in volumes:
        mgr.nxmapset.add_nxmap_handler_from_volume(v)

from chimerax.core.commands.atomspec import AtomSpecArg
class VolumesArg(AtomSpecArg):
    """Parse command models specifier"""
    name = "a models specifier"

    @classmethod
    def parse(cls, text, session):
        '''
        Returns only Volume objects (not subclasses)
        '''
        from chimerax.map import Volume
        aspec, text, rest = super().parse(text, session)
        models = aspec.evaluate(session).models
        volumes = [m for m in models if type(m) == Volume]
        return volumes, text, rest


def register_clipper_cmd(logger):
    from chimerax.core.commands import (
        register, CmdDesc,
        BoolArg, FloatArg
        )
    from chimerax.atomic import StructuresArg, StructureArg, AtomsArg

    spot_desc = CmdDesc(
        optional=[
            ('models', StructuresArg),
            ('enable', BoolArg)
        ],
        synopsis='Switch on/off "Scrolling sphere" visualisation with live atomic symmetry'
    )
    register('cxclipper spotlight', spot_desc, spotlight, logger=logger)

    vol_desc = CmdDesc(
        required=[
            ('volumes', VolumesArg),
        ],
        keyword=[
            ('to_model', StructureArg)
        ],
        synopsis='Have Clipper take control of the chosen volumes and associate them with the given model'
    )
    register('cxclipper associate', vol_desc, associate_volumes, logger=logger)
