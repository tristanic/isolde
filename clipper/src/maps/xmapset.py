from .mapset_base import MapSet_Base

class XmapSet_Live(Mapset_Base):
    '''
    Handles creation, deletion, recalculation and visualisation of
    crystallographic maps based on the current model and a set of observed
    reflection data (FOBS, SIGFOBS).
    '''

    def __init__(self, manager, fsigf_name=None, bsharp_vals=[],
        exclude_free_reflections=False, fill_with_fcalc=False,
        live_update=True, show_r_factors=True):
        '''
        Prepare the C++ :class:`Xtal_mgr` object and create a set of
        crystallographic maps. The standard 2mFo-DFc and mFo-DFc maps will
        always be created, while any number of maps with different
        sharpening/smoothing can be created by providing a list of b-factors in
        bsharp_vals.
        Args:
            manager:
                The master :class:`Map_Mgr`
            fsigf_name:
                The label of the :class:`F_sigF_float` object containing the
                observed amplitudes
            bsharp_vals:
                For each value in this list, a 2mFo-DFc map will be generated
                with the given B_sharp value. A negative B_sharp yields a
                sharpened map, while a positive B_sharp gives smoothing. As a
                rough rule of thumb, a value of B_sharp=-100 tends to be
                reasonable for a 3.8 Angstrom map, and -50 for 3 Angstroms.
                For maps with resolutions better than ~2 Angstroms, it may be
                useful in some circumstances to apply smoothing. This can
                sometimes bring out weaker, low-resolution features hidden by
                high-resolution noise in standard maps.
            exclude_free_reflections:
                If True, observed amplitudes corresponding to the free set will
                not be used in generating the maps. The values used in their
                place will depend on the value of `fill_with_fcalc`. Note that
                this only affects maps generated for viewing - the MDFF potential
                map is always generated with free reflections excluded.
            fill_with_fcalc:
                If `exclude_free_reflections` is False this argument will be
                ignored. Otherwise, if `fill_with_fcalc` is True then the
                excluded amplitudes will be replaced by sigmaa-weighted F(calc).
                If False, the excluded amplitudes will be set to zero.
            display_radius:
                The radius (in Angstroms) of the display sphere used in
                live scrolling mode.
            live_update:
                If True, maps will be automatically recalculated whenever
                coordinates change
            show_r_factors:
                If True, print new R factors to the status bar on each map
                recalculation
        '''
        super().__init__(manager, 'Live real-space maps')
        manager.add([self])

        from chimerax.core.triggerset import TriggerSet
        trig = self._triggers = TriggerSet()

        trigger_names = (
            'maps recalculated'
        )
        for t in trigger_names:
            trig.add_trigger(t)

        self._live_update = False
        self._recalc_needed = False
        self._model_changes_handler = None
        self._delayed_recalc_handler = None
        self._show_r_factors = show_r_factors




    @property
    def triggers(self):
        return self._triggers
