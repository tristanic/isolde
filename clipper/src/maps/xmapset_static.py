from .mapset_base import MapSet_Base
class XmapSet_Static(MapSet_Base):
    '''
    Handles the organisation and visualisation of pre-calculated
    crystallographic maps.
    '''
    def __init__(self, manager, datasets=None):
        super().__init__(manager, 'Pre-calculated crystallographic maps')

    def add_xmap_handler(self, dataset, is_difference_map = None,
                color = None, style = None, contour = None):
        '''
        Add a new XmapHandler based on the given reflections and phases.
        Args:
            dataset:
                a ReflectionData_Calc object.
            is_difference_map:
                Decides whether this map is to be treated as a difference
                map (with positive and negative contours) or a standard
                map with positive contours only. Leave as None to allow
                this to be determined automatically from the dataset,
                or set it to True or False to force an explicit choice.
            color:
                an array of [r, g, b, alpha] as integers in the range 0-255
                for a positive density map, or
                [[r,g,b,alpha],[r,g,b,alpha]] representing the colours of
                the positive and negative contours respectively if the
                map is a difference map
            style:
                one of 'mesh' or 'surface'
            contour:
                The value(s) (in sigma units) at which to contour the map
                on initial loading. For a standard map this should be a
                single value; for a difference map it should be
                [negative contour, positive contour]
        '''
        data = dataset.data
        if is_difference_map is None:
            is_difference_map = dataset.is_difference_map
        if is_difference_map and color is not None and len(color) != 2:
            err_string = '''
            ERROR: For a difference map you need to define colours for
            both positive and negative contours, as:
            [[r,g,b,a],[r,g,b,a]] in order [positive, negative].
            '''
            raise TypeError(err_string)
        new_handler = XmapHandler_Static(self, dataset.name, data,
            *self.box_params,
            is_difference_map = is_difference_map
        )

        if style is None:
            style = 'mesh'
        if color is None:
            if is_difference_map:
                color = self.DEFAULT_DIFF_MAP_COLORS
            elif style == 'mesh':
                color = [self.DEFAULT_MESH_MAP_COLOR]
            else:
                color = [self.DEFAULT_SOLID_MAP_COLOR]
        if contour is None:
            if is_difference_map:
                contour = self.STANDARD_DIFFERENCE_MAP_CONTOURS
            else:
                contour = self.STANDARD_LOW_CONTOUR
        elif not hasattr(contour, '__len__'):
                contour = numpy.array([contour])
        else:
            contour = numpy.array(contour)
        contour = contour * new_handler.sigma
        self.add([new_handler])
        new_handler.set_representation(style)
        new_handler.set_parameters(**{'cap_faces': False,
                                  'surface_levels': contour,
                                  'show_outline_box': False,
                                  'surface_colors': color,
                                  'square_mesh': True})
        return new_handler




from .map_handler_base import XmapHandler_Base
class XmapHandler_Static(XmapHandler_Base):
    '''
    An XmapHandler_Static is in effect a resizable window into a periodic
    crystallographic map. The actual map data (a clipper Xmap object) is
    calculated via fast Fourier transform from a pre-calculated set of
    amplitudes and phases (e.g. as provided by a refinement package), and filled
    into the XmapHandler_Static.data array as needed. Mothods are provided for
    tracking and filling a box around the centre of rotation, and static display
    of a given region.
    '''
    def __init__(self, mapset, name, f_phi_data, origin, grid_origin, grid_dim,
        is_difference_map=False):
        '''
        Args:
            mapset:
                The XmapSet_Live object this belongs to
            name:
                A descriptive name for this map
            f_phi_data:
                The Clipper HKL_Data_F_Phi object the map is to be calculated
                from.
            origin:
                The (x,y,z) coordinates of the bottom left corner of the
                volume.
            grid_origin:
                The (u,v,w) integer grid coordinates corresponding to
                origin.
            grid_dim:
                The shape of the box in (u,v,w) grid coordinates.
            is_difference_map:
                Is this a difference map?
        '''
        super().__init__(mapset, name, origin, grid_origin, grid_dim,
            is_difference_map=is_difference_map)

        from ..clipper_python import Xmap_float as Xmap
        xmap = self._xmap = Xmap(self.spacegroup, self.cell, self.grid_sampling)
        xmap.fft_from(f_phi_data)

    @property
    def xmap(self):
        return self._xmap

    @property
    def stats(self):
        if not hasattr(self, '_stats') or self._stats is None:
            from ..clipper_python import Map_stats
            all = self._all_stats = Map_stats(self)
            self._stats = (all.mean, all.std_dev, all.std_dev)
        return self._stats
