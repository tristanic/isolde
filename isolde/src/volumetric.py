# Copyright 2017 Tristan Croll
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


# vim: set expandtab shiftwidth=4 softtabstop=4:

# TODO: Re-do this to a better standard
class IsoldeMap:
    '''Generic class to hold the details of maps to be used in the simulation'''
    def __init__(self, session, name, source_map, mask_cutoff, 
        coupling_constant, is_difference_map = False, style = None, 
        color = None, contour = None, contour_units = 'sigma', mask = True, 
        crop = True, per_atom_coupling = False, per_atom_scaling_factors = None):
        '''
        Wraps a ChimeraX Volume object (or subclass) with all the added
        information necessary to couple it to an ISOLDE simulation, along
        with some generic methods to adjust its visualisation.
        Args:
            session:
                The ChimeraX session.
            name:
                A descriptive name for this volume to be displayed in the
                GUI etc.
            source_map:
                The existing Volume object
            cutoff:
                The distance (in Angstroms) beyond which the map will be
                masked out during simulations
            coupling_constant:
                A global scaling constant defining how strongly the map 
                pulls on atoms. At each timepoint, atom i feels a force
                equal to:
                    
                    per_atom_scaling_factors[i] * coupling_constant * (map gradient)
                
                along each axis.
            is_difference_map:
                Is this map a difference map? If so, it will be displayed
                with symmetric positive and negative contours.
            style:
                Contour display style for this map. One of "mesh" or "surface"
            color:
                Colour of the map surface. For a standard map, this should
                be a single array [r,g,b,opacity], with each entry in the
                range 0-255. For a difference map, a tuple of two such 
                arrays is required, with the colour for the negative 
                contour first.
            contour:
                The desired contour level of the map display, in units 
                defined by the contour_units argument. Default is 1.0 
                sigma.
            contour_units:
                'sigma' (default): contour will be interpreted as a 
                    multiple of the map standard deviation
                'map units': contour will be interpreted as an absolute
                    value in whatever units the map uses.
            mask:
                boolean (default True). Do we want to mask the map to 
                the mobile atoms during a simulation?
            crop:
                boolean (default True). Do we need to crop down from a 
                bigger map before starting the simulation? If true, the
                default padding will be set to twice the mask radius.
            per_atom_coupling:
                (NOT YET IMPLEMENTED)
                boolean (default False). If true, then the coupling of 
                each atom to the map will be scaled according to its entry
                in per_atom_scaling_factors
            
                
        '''
        self.session = session # Handle for current ChimeraX session
        self._name = name     # User-specified name (e.g. '2mFo-DFc')
        self._source_map = source_map # A model currently loaded into ChimeraX
        self._mask_cutoff = mask_cutoff # in Angstroms 
        self._coupling_constant = coupling_constant # How hard map pulls on atoms
        self._is_difference_map = is_difference_map # Is this a Fo-Fc or similar difference map?
        self._per_atom_coupling = per_atom_coupling # Do we vary map pull by atom type?
        self._style = style # Map surface style
        self._color = color # Map display color
        self._contour = contour # Map contour level
        self._contour_units = contour_units # Units for contour level (sigma or map units)
        self._mask = mask # Do we mask down the volume to the mobile selection?
        # TODO: Add the ability to define specific sub-selections of atoms
        # to be associated with a particular map (e.g. anomalous scatterers with
        # an anomalous difference map; omitted fragments with the mFo-DFc density
        # in an omit map, etc.)
        self._crop = crop # Do we crop the map before starting the simulation?
        
        # Optionally, we can scale the strength of each atom's individual coupling to the map
        if per_atom_coupling:
            # dict relating atom index to scale factor
            self._per_atom_coupling_scale_factor = {}
        else:
            # single global value
            self._per_atom_coupling_scale_factor = 1.0
        
        
        self._source_map_res = source_map.region_origin_and_step(source_map.region)[1]
        # TODO: This currently ignores the skewness of the map (in most cases
        # probably not a huge deal). Still, it should be re-calculated to give
        # the effective resolution along the model (x,y,z) axes
        
        # Placeholder for the masked map
        self._masked_map = None        
        # Placeholder for the Continuous3DFunction that will be generated by OpenMM
        self._c3d_function = None
        
        # Placeholder for the final potential force function generated by OpenMM
        self._potential_function = None
        
    
    def change_map_parameters(self, source_map = None, 
                                cutoff = None, 
                                coupling_constant = None,
                                is_difference_map = None,
                                style = None,
                                color = None,
                                contour = None,
                                contour_units = None,
                                mask = None,
                                crop = None,
                                per_atom_coupling = None 
                                ):
        if source_map is not None: 
            self._source_map = source_map 
        if cutoff is not None: 
            self._mask_cutoff = cutoff 
        if coupling_constant is not None: 
            self._coupling_constant = coupling_constant 
        if style is not None: 
            self._style = style
        if color is not None: 
            self._color = color 
        if contour is not None: self._contour = contour 
        if contour_units is not None: 
            self._contour_units = contour_units 
        if mask is not None: 
            self._mask = mask 
        if crop is not None:
            self._crop = crop
        if per_atom_coupling is not None: 
            self._per_atom_coupling = per_atom_coupling 
        
    
    def set_source_map(self, source_map):
        self._source_map = source_map
    
    # Boolean switch: do we mask to the mobile atoms when the simulation starts?
    def set_mask_vis(self, mask):
        self._mask = mask

    @property
    def mask_cutoff(self):
        ''' 
        Cutoff in Angstroms beyond which the map will be masked out 
        during simulations.
        '''
        return self._mask_cutoff
    
    @mask_cutoff.setter
    def mask_cutoff(self, cutoff):
        self._mask_cutoff = cutoff
    
    def set_mask_cutoff(self, cutoff):
        self._mask_cutoff = cutoff
        
    def set_display_style(self, style):
        if style not in ['mesh', 'surface', 'solid']:
            print('style must be one of \'mesh\', \'surface\' or \'solid\'')
            return
        else:
            self._style = style
    
    def set_color(self, color1_name, color2_name = None):
        if self._is_difference_map:
            if color2_name is None:
                raise TypeError('A difference map requires two colours!')
        else:
            if color2_name is not None:
                raise TypeError('A standard map has only one colour!')
        from chimerax.core.map import volumecommand
        from chimerax.core import colors
        rgba = [colors.BuiltinColors[color1_name]]
        if color2_name is not None:
            rgba.append(colors.BuiltinColors[color2_name])
        volumecommand.volume(session, [self._source_map], color = [rgba])
        self._color = [color1_name, color2_name]
        
    def set_contour(self, contour):
        self._contour = contour
        
    def set_contour_units(self, text):
        if text not in ['sigma', 'map units']:
            print('text must be either \'sigma\' or \'map units\'')
            return
        else:
            self._contour_units = text
    
    
    def set_coupling_constant(self, coupling_constant):
        self._coupling_constant = coupling_constant
        
    def set_c3d_function(self, func):
        self._c3d_function = func
        
    def set_potential_function(self, func):
        self._potential_function = func
    
    def get_contour(self):
        return self._contour, self._contour_units
    
    def get_map_parameters(self):
        return self._name, self._source_map, self._mask_cutoff, \
                self._coupling_constant, self._style, self._color, \
                self._contour, self._contour_units, self._mask, \
                self._per_atom_coupling, self._per_atom_coupling_scale_factor
    
    def per_atom_coupling(self):
        return self._per_atom_coupling
    
    def get_per_atom_coupling_params(self):
        return self._per_atom_coupling_scale_factor
    
    
    # Set the per-atom scale factor for a single atom    
    def set_per_atom_coupling_constant(self, index, value):
        self._per_atom_coupling_scale_factor[index] = value
        
    # Set all per-atom scale factors at once (e.g. to copy parameters from another map)
    def set_per_atom_coupling_constants(self, params):
        from copy import deepcopy
        self._per_atom_coupling_scale_factor = deepcopy(params)
                
    def get_name(self):
        return self._name
    
    def get_source_map(self):
        return self._source_map

    def get_mask_vis(self):
        return self._mask
    
    def get_mask_cutoff(self):
        return self._mask_cutoff
    
    def get_color(self):
        return self._color
    
    def get_display_style(self):
        return self._style
    
    def get_coupling_constant(self):
        return self._coupling_constant
        
    def get_c3d_function(self):
        return self._c3d_function
        
    def get_potential_function(self):
        return self._potential_function
        
    coupling_constant = property(get_coupling_constant, set_coupling_constant)
    
    def crop_to_selection(self, selection, padding, normalize = False):
        '''
        Crop a map to surround a selection, plus the given padding (in 
        Angstroms). 
        Args:
            selection:
                an Atoms object
            padding:
                the minimum distance of the box edge from any atom
            normalize:
                if true, all values within the returned map will be divided
                by the overall standard deviation of the parent map.
        '''
        import numpy
        if self._crop:
            points = selection.scene_coords
            self._source_map.position.inverse().move(points) # Convert points to map coordinates
            from chimerax.core.map.data.regions import points_ijk_bounds, \
                clamp_region, integer_region
            data = self._source_map.data
            r = points_ijk_bounds(points, padding, data)
            r = clamp_region(integer_region(r), data.size)
            from chimerax.core.map.data.griddata import Grid_Subregion
            grid_data = Grid_Subregion(data, r[0],r[1])
            m = numpy.empty(grid_data.matrix().shape, numpy.float64)
            m[:] = grid_data.matrix()
        else:
            grid_data = self._source_map.data
            orig_data = grid_data.matrix()
            m = numpy.empty(orig_data.shape, numpy.float64)
            m[:] = orig_data
        if normalize:
            m /= self._source_map.mean_sd_rms()[1]
        return m, grid_data
    
    def mask_volume_to_selection(self, selection, resolution = None, invert = False, normalize = False):
        '''
        Mask a given map down to a selection and interpolate it onto an orthogonal
        grid of the specified resolution. Resolution must be either a single
        number or an (x,y,z) numpy array. If resolution is not specified, 
        the resolution of the source map will be used. Optionally, the map 
        may also be inverted or normalised such that its standard deviation = 1.
        '''

        big_map = self._source_map
        cutoff = self._mask_cutoff
        sel = selection
        import numpy as np
        # Get minimum and maximum coordinates of selected atoms, and pad
        maxcoor = (sel.coords.max(0)  + cutoff)
        mincoor = (sel.coords.min(0)  - cutoff)
        # ChimeraX wants the volume array to be in zyx order, so we need to reverse
        vol_size = maxcoor[::-1] - mincoor[::-1]
        # Round to an integral number of steps at the desired resolution
        if resolution is not None:
            res = resolution
        else:
            res = big_map.data_origin_and_step()[1]
        sizek, sizej, sizei = np.ceil(vol_size/res).astype(int)
        vol_array = np.zeros([sizek, sizej, sizei])
        
        from chimerax.core.map.data import Array_Grid_Data
        from chimerax.core.map import Volume
        
        mask_array = Array_Grid_Data(vol_array, origin = mincoor, step = res*np.ones(3))
        vol = Volume(mask_array, self.session)
        vol_coords = vol.grid_points(vol.model_transform())
        
        from chimerax.core.geometry import find_close_points
        map_points, ignore = find_close_points(vol_coords, sel.coords, cutoff)
        mask_1d = np.zeros(len(vol_coords))
        mask_1d[map_points] = 1
        vol_array = np.reshape(mask_1d, (sizek, sizej, sizei))
        if normalize:
            std = big_map.mean_sd_rms()[1]
            cropped_map = big_map.interpolate_on_grid(vol) / std
        else:
            cropped_map = big_map.interpolate_on_grid(vol)
        if invert:
            masked_map = -cropped_map[0] * vol_array
        else:
            masked_map = cropped_map[0] * vol_array
        masked_array = Array_Grid_Data(masked_map, origin = mincoor, step = res*np.ones(3))
        vol = self._masked_map = Volume(masked_array, self.session)
        return vol

