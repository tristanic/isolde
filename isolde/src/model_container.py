# Top-level container for one ISOLDE modelling target

from chimerax.core.models import Model

class Isolde_Model_Container(Model):
    '''
    A model prepared for use in ISOLDE typically has a great deal of data
    associated with it: the atomic coordinates themselves, crystallographic
    datasets and map(s), real-space map(s), reference models, etc.. The design
    of various core features of ChimeraX dictates that many such objects (in
    particular, atomic structures and non-crystallographic maps) must not be
    subordinate to the main atomic model. Therefore, we need a top-level model
    to act as a container for data organisation.
    '''
    def __init__(self, session, isolde, structure,
            reference_models=[],
            real_space_maps=[],
            precalculated_map_mtz=None,
            xtal_data_mtz=None, live_xtal_map_recalc=True):
        name = "ISOLDE: "+structure.name
        super().__init__(name, session)
        self.isolde = isolde
        self.add([structure])
        self._atomic_model = structure

        if len(reference_models):
            raise NotImplementedError('Reference models not yet supported')
        if len(real_space_maps):
            raise NotImplementedError('Real-space maps not yet supported')
        if precalculated_map_mtz:
            raise NotImplementedError('Precalculated crystallographic maps not yet supported')

        if xtal_data_mtz is not None:
            from .isolde import initialize_live_xtal_structure
            initialize_live_xtal_structure(xtal_data_mtz, live=live_xtal_map_recalc)

        session.models.add([self])

    @property
    def atomic_model(self):
        return self._atomic_model

    def initialize_live_xtal_structure(self, filename, live=True):
        m = self.atomic_model
        from chimerax.clipper import symmetry
        sh = symmetry.XtalSymmetryHandler(m, mtzfile=filename,
            map_oversampling=self.params.map_shannon_rate)
        sh.xmapset.live_update = live
        standard_map = sh.xmapset['2mFo-DFc']
        diff_map = sh.xmapset['mFo-DFc']
        has_extra_map = False
        extra_map_is_sharp = False
        extra_map = None
        for key, xmap in sh.xmapset.items():
            if "sharp" in key:
                has_extra_map = True
                extra_map_is_sharp = True
                extra_map = xmap
            elif "smooth" in key:
                has_extra_map = True
                extra_map = xmap
        from . import visualisation as v
        from chimerax.map import volumecommand
        sd = standard_map.mean_sd_rms()[1]
        diff_styleargs = v.map_style_settings[v.map_styles.solid_t40]
        volumecommand.volume(self.session, [diff_map], **diff_styleargs)
        if has_extra_map:
            xsd = extra_map.mean_sd_rms()[1]
            styleargs = v.map_style_settings[v.map_styles.solid_t20]
            if extra_map_is_sharp:
                volumecommand.volume(self.session, [extra_map], **styleargs)
                standard_map.set_parameters(surface_levels = (1.5*sd,))
                extra_map.set_parameters(surface_levels=(3.0*xsd,))
            else:
                volumecommand.volume(self.session, [standard_map], **styleargs)
                standard_map.set_parameters(surface_levels = (3.0*sd,))
                extra_map.set_parameters(surface_levels=(1.5*xsd,))
        # Set difference map according to standard map SD, since this is more
        # meaningful
        diff_map.set_parameters(surface_levels=(-0.8*sd, 0.8*sd))

        # Set the MDFF coupling constant
        mdff_p = sh.xmapset['MDFF potential']
        from .session_extensions import get_mdff_mgr
        mdff_mgr = get_mdff_mgr(m, mdff_p)
        from .openmm.weighting import guess_mdff_weight
        sim_params = self.isolde.sim_params
        mdff_mgr.global_k = guess_mdff_weight(
            mdff_mgr,
            constant=sim_params.standard_map_coupling_base_constant
        )

    @property
    def mdff_potentials(self):
        '''
        Returns a list of 2-tuples of :class:`Volume` (or subclass) and
        :class:`MDFF_Mgr` instances containing all volumes currently acting as
        MDFF potentials.
        '''
        potentials = []
        from .session_extensions import get_mdff_mgr
        from chimerax.map import Volume
        for m in self.all_models:
            if isinstance(m, Volume):
                mgr = get_mdff_mgr(self.atomic_model, m)
                if mgr is not None:
                    potentials.append((m, mgr))
        return potentials
