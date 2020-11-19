# @Author: Tristan Croll <tic20>
# @Date:   29-Oct-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 02-Nov-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

class Water_Finder:
    def __init__(self, session, structure, surface):
        self.session=session
        self.structure = structure
        self.surface=surface
        self._blobs = []
        self._blob_centers = None

    def find_blobs(self, enclosed_min_volume=0.1, enclosed_max_volume=8.0):
        self._blobs = []
        #v = self.volume
        surface = self.surface
        # if surface is None:
        #     if len(v.surfaces==0):
        #         raise RuntimeError('Volume has no surface defined! Add a surface '
        #             'with a contour level showing waters first.')
        #     elif len(v.surfaces > 1):
        #         raise RuntimeError('Volume has multiple surfaces defined! Please '
        #             'specify explicitly which one you wish to use for water '
        #             'picking.')
        #     surface = v.surfaces
        from chimerax.surface import connected_pieces
        from chimerax.surface.area import enclosed_volume
        triangles = surface.masked_triangles
        vertices = surface.vertices
        blob_list = connected_pieces(triangles)
        for (vi, ti) in blob_list:
            vol, holes = enclosed_volume(vertices, triangles[ti,:])
            if vol is None or vol > enclosed_max_volume or vol < enclosed_min_volume:
                continue
            self._blobs.append((vi, ti))

    def blob_centers(self):
        surface = self.surface
        vertices = surface.vertices
        from chimerax.surface import vertex_areas
        triangles = surface.masked_triangles
        blob_list = self._blobs

        varea = vertex_areas(vertices, triangles)
        from numpy import empty, float32
        centers = empty((len(blob_list), 3), float32)
        for i, (vi, ti) in enumerate(blob_list):
            blob_varea = varea[vi]
            blob_area = blob_varea.sum()
            centers[i] = blob_varea.dot(vertices[vi])/blob_area
        self._blob_centers = centers

    def find_plausible_waters(self, min_distance=2.2, max_distance=3.2):
        # TODO: base this on actual van der Waals radii
        m = self.structure
        # For now, let's treat anything other than carbon or hydrogen as a
        # potential h-bonder
        nonh = m.atoms[m.atoms.element_names != 'H']
        h_bonders = nonh[nonh.element_names != 'C']
        from chimerax.geometry import find_close_points
        bc = self._blob_centers
        all_coords = nonh.coords
        a_coords = h_bonders.coords
        too_close, _ = find_close_points(bc, all_coords, min_distance)
        import numpy
        bc = numpy.delete(bc, too_close, axis=0)
        if len(bc):
            water_indices, _ = find_close_points(bc, a_coords, max_distance)
            bc = bc[water_indices]
        if not len(bc):
            return
        from chimerax.isolde.atomic.building.place_ligand import place_water
        new_waters = []
        for pos in bc:
            w = place_water(self.session, m, pos, distance_cutoff=max_distance,
                sim_settle=False)
            new_waters.append(w)
        return new_waters
