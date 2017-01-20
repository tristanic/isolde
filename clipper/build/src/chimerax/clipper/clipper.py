from .lib import clipper_python_core as clipper_core
import numpy

class Atom():
    '''
    A minimalist atom object containing only the properties required for
    electron density calculations
    '''
    ATOM_NAMES = ['H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',
                  'Ne', 'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar',
                  'K',  'Ca', 'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co',
                  'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
                  'Rb', 'Sr', 'Y',  'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',
                  'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',  'Xe',
                  'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu',
                  'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
                  'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl',
                  'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
                  'Pa', 'U',  'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 
                  'H1-',  'Li1+', 'Be2+', 'Cval', 'O1-',  'O2-',  'F1-',
                  'Na1+', 'Mg2+', 'Al3+', 'Siva', 'Si4+', 'Cl1-', 'K1+',
                  'Ca2+', 'Sc3+', 'Ti2+', 'Ti3+', 'Ti4+', 'V2+',  'V3+',
                  'V5+',  'Cr2+', 'Cr3+', 'Mn2+', 'Mn3+', 'Mn4+', 'Fe2+',
                  'Fe3+', 'Co2+', 'Co3+', 'Ni2+', 'Ni3+', 'Cu1+', 'Cu2+',
                  'Zn2+', 'Ga3+', 'Ge4+', 'Br1-', 'Rb1+', 'Sr2+', 'Y3+',
                  'Zr4+', 'Nb3+', 'Nb5+', 'Mo3+', 'Mo5+', 'Mo6+', 'Ru3+',
                  'Ru4+', 'Rh3+', 'Rh4+', 'Pd2+', 'Pd4+', 'Ag1+', 'Ag2+',
                  'Cd2+', 'In3+', 'Sn2+', 'Sn4+', 'Sb3+', 'Sb5+', 'I1-',
                  'Cs1+', 'Ba2+', 'La3+', 'Ce3+', 'Ce4+', 'Pr3+', 'Pr4+',
                  'Nd3+', 'Pm3+', 'Sm3+', 'Eu2+', 'Eu3+', 'Gd3+', 'Tb3+',
                  'Dy3+', 'Ho3+', 'Er3+', 'Tm3+', 'Yb2+', 'Yb3+', 'Lu3+',
                  'Hf4+', 'Ta5+', 'W6+',  'Os4+', 'Ir3+', 'Ir4+', 'Pt2+',
                  'Pt4+', 'Au1+', 'Au3+', 'Hg1+', 'Hg2+', 'Tl1+', 'Tl3+',
                  'Pb2+', 'Pb4+', 'Bi3+', 'Bi5+', 'Ra2+', 'Ac3+', 'Th4+',
                  'U3+',  'U4+',  'U6+',  'Np3+', 'Np4+', 'Np6+', 'Pu3+',
                  'Pu4+', 'Pu6+']  
                        
    def __init__(self):
        self._core_atom = clipper_core.Atom()
    
    @property
    def core_atom(self):
        return self._core_atom
    
    @property
    def element(self):
        '''
        The standard abbreviated element name. All valid names are listed in
        Atom.ATOM_NAMES
        '''
        return self.core_atom.element()
    
    @element.setter
    def element(self, element_name):
        # Check common atom names first to improve performance
        if element_name not in ('H', 'C', 'N', 'O', 'S'):
            if element_name not in self.ATOM_NAMES:
                raise TypeError('Unrecognised element!')
        self.core_atom.set_element(element_name)
    
    @property
    def coord(self):
        '''
        Get the coordinates of the atom as a numpy array
        '''
        coord_orth = self.core_atom.coord_orth()
        return numpy.array([coord_orth.x(), coord_orth.y(), coord_orth.z()])
    
    @coord.setter
    def coord(self, coord):
        '''
        Set the coordinates of the atom using a list or numpy array
        '''
        self.core_atom.set_coord_orth(clipper_core.Coord_orth(*coord))

    @property
    def coord_orth(self):
        '''
        Get the Clipper Coord_orth object associated with this atom
        '''
        return self.core_atom.coord_orth()
    
    @coord_orth.setter
    def coord_orth(self, coord):
        self.coord = coord
        
    @property
    def occupancy(self):
        return self.core_atom.occupancy()
    
    @occupancy.setter
    def occupancy(self, occ):
        self.core_atom.set_occupancy(occ)
    
    @property
    def u_iso(self):
        '''
        Get the isotropic b-factor
        '''
        return self.core_atom.u_iso()
    
    @u_iso.setter
    def u_iso(self, u_iso):
        self.core_atom.set_u_iso(u_iso)
    
    @property
    def b_factor(self):
        '''
        Synonym for u_iso
        '''
        return self.u_iso
    
    @b_factor.setter
    def b_factor(self, b_factor):
        self.u_iso = b_factor
    
    @property
    def u_aniso(self):
        '''
        Get the object holding the anisotropic B-factor matrix. Note that
        Clipper-python does not currently have functions to retrieve the
        matrix elements
        '''
        return self.core_atom.u_aniso_orth()
    
    @u_aniso.setter
    def u_aniso(self, u11, u22, u33, u12, u13, u23):
        self.core_atom.set_u_aniso_orth(u11, u22, u33, u12, u13, u23)
    
    @property
    def is_null(self):
        return self.core_atom.is_null()
    
    def transform(self, rot_trans):
        '''
        Transform the atom using the rotation and translation defined in
        a Clipper RTop_orth object
        '''
        self.core_atom.transform(rot_rans)

class Xmap(clipper_core.Xmap_double):
    def __init__(self, session, name, spacegroup, cell, grid_sam):
        clipper_core.Xmap_double.__init__(self, spacegroup, cell, grid_sam)
        self.session = session
        self.name = name
        # Some extra useful variables that aren't directly available from
        # the Clipper API
        c = self.cell()
        g = self.grid_sampling()
        self.grid_samples = g.dim()
        
        ###
        # Basic stats. Can only be set after the map has been filled with an
        # FFT. Returned as a tuple in the below order by self.stats()
        ###
        self.max, self.min, self.mean, \
            self.sigma, self.skewness, self.kurtosis = self.stats()
        
        # Default "radius" (actually half-width of the rhombohedron) of
        # the box in which the map will be drawn.
        self.box_radius_angstroms = 20
        # Box dimensions in grid points (depends on resolution)
        self.box_dimensions = None
        # Centre point of the box (will typically be set to the centre of
        # rotation).
        self.box_center = None
        # Vector mapping the box centre to its origin
        self.box_origin_offset = None
        # ChimeraX Volume object to draw the map into
        self.volume = None
        # Array_Grid_Data object held by the Volume object
        self.array_grid_data = None
        # Numpy array to send the map data to
        self.box_data = None
        # Distance any axis of the cofr has to move before triggering an update
        self.cofr_eps = 0.01
        
        # Model object to hold Drawings defining the special positions
        self.special_positions_model = None
        
    def calculate_stats(self):
        self.max, self.min, self.mean, \
            self.sigma, self.skewness, self.kurtosis = self.stats()         
    
    def initialize_box_display(self, radius = 15, pad = 1):
        from chimerax.core.commands import camera, cofr
        camera.camera(self.session, 'ortho')
        cofr.cofr(self.session, 'centerOfView')
        self.box_radius = radius
        self.box_pad = pad
        self.box_radius_angstroms = radius
        v = self.session.view
        c = self.cell()
        g = self.grid_sampling()
        self.box_center = v.center_of_rotation
        box_corner_grid, box_corner_xyz = self.find_box_corner(self.box_center, radius, pad)
        self.box_origin_offset = box_corner_xyz - self.box_center
        # Add padding of one voxel on each side to make sure we're always above
        # the radius
        self.box_dimensions = (numpy.ceil(radius / self.voxel_size() * 2)+pad*2).astype(int)[::-1]
        data = self.box_array = numpy.ones(self.box_dimensions, numpy.double)
                
        
        from chimerax.core.map.data import Array_Grid_Data
        from chimerax.core.map import Volume
        self.array_grid_data = Array_Grid_Data(data, origin = box_corner_xyz,
                step = self.voxel_size(), cell_angles = c.angles_deg())
        self.volume = Volume(self.array_grid_data, self.session)
        self.fill_volume_data(data, box_corner_grid)
        self.volume.initialize_thresholds()
        self.session.models.add([self.volume])
        contour_val = self.contour = 1.5 * self.sigma
        self.change_contour(contour_val)
        self.volume.set_color([0.5,1.0,0.5,0.6])
        self.volume.show()
        self.session.triggers.add_handler('new frame', self.update_box)
        
    def update_box(self, *_):
        v = self.session.view
        cofr = v.center_of_rotation
        if numpy.all(abs(self.box_center - cofr) < self.cofr_eps):
            return
        self.box_center = cofr
        box_corner_grid, box_corner_xyz = self.find_box_corner(cofr, self.box_radius, self.box_pad)
        self.array_grid_data.set_origin(box_corner_xyz)
        self.fill_volume_data(self.box_array, box_corner_grid)
        self.change_contour(self.contour)
        
    def change_contour(self, contour_val):
        from chimerax.core.map import volumecommand
        volumecommand.volume(self.session,[self.volume], level=[[contour_val]], cap_faces = False)
        

    def fill_volume_data(self, target, start_grid_coor):
        shape = (numpy.array(target.shape)[::-1] - 1).tolist()
        end_grid_coor = start_grid_coor + clipper_core.Coord_grid(*shape)
        count = self.export_section_numpy(target, start_grid_coor, end_grid_coor, 'C')
                
    
    def find_box_corner(self, center, radius = 20, pad = 1):
        '''
        Find the bottom corner (i.e. the origin) of a rhombohedral box
        big enough to hold a sphere of the desired radius, and padded by
        pad voxels in each dimension.
        '''
        cell = self.cell()
        grid = self.grid_sampling()
        voxel_size = self.voxel_size()
        radii_in_voxels = radius / voxel_size
        radii_frac = clipper_core.Coord_frac(*(radii_in_voxels * self.voxel_size_frac()))
        center_ortho = clipper_core.Coord_orth(*center)
        center_frac = center_ortho.coord_frac(cell)
        bottom_corner_frac = center_frac - radii_frac
        bottom_corner_grid = bottom_corner_frac.coord_grid(grid) -\
                            clipper_core.Coord_grid(pad,pad,pad)
        bottom_corner_orth = bottom_corner_grid.coord_frac(grid).coord_orth(cell)
        return bottom_corner_grid, bottom_corner_orth.xyz()
        
    def draw_special_positions(self):
        from chimerax.core.models import Model, Drawing
        from chimerax.core.geometry import Place
        from chimerax.core.surface.shapes import sphere_geometry
        import copy
        
        m = self.special_positions_model
        
        if m is None or m.deleted:
            m = self.special_positions_model = Model('Special Positions',self.session)
            spc = numpy.array(self.special_positions_unit_cell(self.grid_sampling()))
            sphere = numpy.array(sphere_geometry(4))
            sphere[0]*=0.5
            
            for row in spc:
                coords = row[0:3].tolist()
                mult = row[3]
                cg = clipper_core.Coord_grid(*coords)
                cf = cg.coord_frac(self.grid_sampling())
                co = cf.coord_orth(self.cell())
                d = Drawing('point')
                d.vertices, d.normals, d.triangles = copy.copy(sphere)
                if mult == 3:
                    d.name = '3-fold'
                    d.vertices = d.vertices * 1.5
                    d.set_color([0,255,255,255])
                elif mult == 4:
                    d.name = '4-fold'
                    d.vertices = d.vertices * 2
                    d.set_color([255,255,0,255])
                elif mult == 6:
                    d.name = '6-fold'
                    d.vertices = d.vertices * 3
                    d.set_color([255,0,0,255])
                else:
                    d.name = '2-fold'
                d.set_position(Place(origin = co.xyz()))
                m.add_drawing(d) 
            self.session.models.add([m])
                            
    def draw_special_positions(self):
        from chimerax.core.models import Model, Drawing
        from chimerax.core.geometry import Place, Places
        from chimerax.core.surface.shapes import sphere_geometry
        import copy
        
        m = self.special_positions_model
        
        if m is None or m.deleted:
            m = self.special_positions_model = Model('Special Positions',self.session)
            spc = numpy.array(self.special_positions_unit_cell(self.grid_sampling()))
            sphere = numpy.array(sphere_geometry(80))
            sphere[0]*=0.5
            scale_2fold = numpy.identity(3)
            scale_3fold = numpy.identity(3)* 1.5
            scale_4fold = numpy.identity(3)* 2
            scale_6fold = numpy.identity(3)* 3
            rgba_2fold = numpy.array([255,255,255,255],numpy.int32)
            rgba_3fold = numpy.array([0,255,255,255],numpy.int32)
            rgba_4fold = numpy.array([255,255,0,255],numpy.int32)
            rgba_6fold = numpy.array([255,0,0,255],numpy.int32)
            positions = []
            colors = []
            d = Drawing('points')
            d.vertices, d.normals, d.triangles = sphere
            
            for row in spc:
                coords = row[0:3].tolist()
                mult = row[3]
                cg = clipper_core.Coord_grid(*coords)
                cf = cg.coord_frac(self.grid_sampling())
                co = cf.coord_orth(self.cell())
                if mult == 2:
                    positions.append(Place(axes=scale_2fold, origin=co.xyz()))
                    colors.append(rgba_2fold)
                if mult == 3:
                    positions.append(Place(axes=scale_3fold, origin=co.xyz()))
                    colors.append(rgba_3fold)
                elif mult == 4:
                    positions.append(Place(axes=scale_4fold, origin=co.xyz()))
                    colors.append(rgba_4fold)
                elif mult == 6:
                    positions.append(Place(axes=scale_6fold, origin=co.xyz()))
                    colors.append(rgba_6fold)
            d.positions = Places(positions)
            d.colors = numpy.array(colors)
            m.add_drawing(d)
            self.session.models.add([m])
        
    def draw_special_positions(self, model, offset = None):
        from chimerax.core.models import Model, Drawing
        from chimerax.core.geometry import Place, Places
        from chimerax.core.surface.shapes import sphere_geometry
        import copy
        
        ref = model.bounds().center().astype(float)
        frac_coords = clipper_core.Coord_orth(*ref).coord_frac(self.cell()).uvw()
        if offset is None:
            offset = numpy.array([0,0,0],int)
        corners_frac = numpy.array([[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]],numpy.double) + offset
        corners = []
        for c in corners_frac:
            cf = clipper_core.Coord_frac(*c).coord_orth(self.cell())
            corners.append(cf.xyz())
        m = self.special_positions_model
        
        if m is None or m.deleted:
            m = self.special_positions_model = Model('Special Positions',self.session)
            spc = numpy.array(self.special_positions_unit_cell_xyz(offset))
            coords = spc[:,0:3]
            multiplicity = spc[:,3].astype(int)
            sphere = numpy.array(sphere_geometry(80))
            sphere[0]*=0.5
            scale_2fold = numpy.identity(3)
            scale_3fold = numpy.identity(3)* 1.5
            scale_4fold = numpy.identity(3)* 2
            scale_6fold = numpy.identity(3)* 3
            rgba_2fold = numpy.array([255,255,255,255],numpy.int32)
            rgba_3fold = numpy.array([0,255,255,255],numpy.int32)
            rgba_4fold = numpy.array([255,255,0,255],numpy.int32)
            rgba_6fold = numpy.array([255,0,0,255],numpy.int32)
            rgba_corner = numpy.array([255,0,255,255],numpy.int32)
            positions = []
            colors = []
            d = Drawing('points')
            d.vertices, d.normals, d.triangles = sphere
            
            for coord, mult in zip(coords, multiplicity):
                if mult == 2:
                    positions.append(Place(axes=scale_2fold, origin=coord))
                    colors.append(rgba_2fold)
                elif mult == 3:
                    positions.append(Place(axes=scale_3fold, origin=coord))
                    colors.append(rgba_3fold)
                elif mult == 4:
                    positions.append(Place(axes=scale_4fold, origin=coord))
                    colors.append(rgba_4fold)
                elif mult == 6:
                    positions.append(Place(axes=scale_6fold, origin=coord))
                    colors.append(rgba_6fold)
            for c in corners:
                positions.append(Place(axes=scale_6fold, origin=c))
                colors.append(rgba_corner)
                
            d.positions = Places(positions)
            d.colors = numpy.array(colors)
            m.add_drawing(d)
            self.session.models.add([m])
        
         
class HKL_info(clipper_core.HKL_info):
    def __init__(self, session):
        clipper_core.HKL_info.__init__(self)
        self.session = session
                
        
        
    
    



def apply_b_factors_to_hydrogens(atom_list):
    '''
    When hydrogens are added to a pre-existing structure, their isotropic
    and/or anisotropic B-factors are not necessarily set. This routine simply
    sets them to match their adjacent heavy atoms. The input atom_list can
    be either the hydrogens alone, or the whole structure.
    '''
    hydrogens = atom_list.filter(atom_list.element_names == 'H')
    # Get heavy atoms by taking advantage of the fact that hydrogens are
    # only ever covalently bonded to a single other atom.
    for h in hydrogens:
        b = h.bonds[0]
        for atom in b.atoms:
            if atom.element_name != 'H':
                h.bfactor = atom.bfactor
                h.aniso_u6 = atom.aniso_u6
                break
                
def Atom_list_from_sel(atom_list):
    n = len(atom_list)
    elements = atom_list.element_names
    coords = atom_list.coords
    occupancies = atom_list.occupancy.astype(float)
    u_iso = atom_list.bfactors.astype(float)
    #u_aniso = sel.aniso_u
    u_aniso = [numpy.random.rand(3,3).astype(float)]*n # FIXME For testing only
    
    clipper_atom_list = clipper_core.Atom_list(n)
    
    for i in range(n):
        thisatom = clipper_atom_list[i]
        thisatom.set_element(elements[i])
        thisatom.set_coord_orth(clipper_core.Coord_orth(*coords[i]))
        thisatom.set_occupancy(occupancies[i])
        thisatom.set_u_iso(u_iso[i])
        ua = u_aniso[i]
        thisatom.set_u_aniso_orth(
            clipper_core.U_aniso_orth(
                ua[0][0], ua[1][1], ua[2][2], ua[0][1], ua[0][2], ua[1][2]))
    
    return clipper_atom_list
 
def Atom_list_from_sel_fast(atom_list):
    '''
    Takes a ChimeraX Atoms object, and creates a Clipper Atoms_list object
    from the relevant atom properties.
    '''
    n = len(atom_list)
    elements = atom_list.element_names.tolist()
    coords = atom_list.coords
    occupancies = atom_list.occupancy
    u_iso = atom_list.bfactors
    # Have to be careful here. Atoms.aniso_u6 returns None if any atom in
    # the array has no aniso_u6 entry. 
    u_aniso = atom_list.aniso_u6
    if u_aniso is None:
        u_aniso = numpy.zeros([n,6],numpy.float32)
        # FIXME Once new ChimeraX build arrives with Atoms.has_aniso_u entry
        for i, a in enumerate(atom_list):
            if a.aniso_u6 is not None:
                u_aniso[i] = a.aniso_u6
    clipper_atom_list = clipper_core.Atom_list()
    clipper_atom_list.extend_by(n)

    clipper_atom_list.set_elements(elements)
    clipper_atom_list.set_coord_orth(coords)
    clipper_atom_list.set_occupancies(occupancies)
    clipper_atom_list.set_u_isos(u_iso)
    clipper_atom_list.set_u_anisos(u_aniso)
    return clipper_atom_list

def import_Xmap_from_mtz_test(session, filename):
    myhkl = HKL_info(session)
    fphidata =  clipper_core.HKL_data_F_phi_double()  
    mtzin = clipper_core.CCP4MTZfile()
    mtzin.open_read(filename)
    mtzin.import_hkl_info(myhkl)
    mtzin.import_hkl_data(fphidata, '/crystal/dataset/[2FOFCWT, PH2FOFCWT]')
    mtzin.close_read()
    name = '2FOFCWT'
    mygrid = clipper_core.Grid_sampling(myhkl.spacegroup(), myhkl.cell(), myhkl.resolution())
    #mymap = clipper_core.Xmap_double(myhkl.spacegroup(), myhkl.cell(), mygrid)
    mymap = Xmap(session, name, myhkl.spacegroup(), myhkl.cell(), mygrid)
    mymap.fft_from(fphidata)
    mymap.calculate_stats()
    return (myhkl, mymap)
    
def vol_box(hklinfo, xmap, min_coor, max_coor):
    cell = hklifno.cell()
    grid = xmap.grid_sampling()
    min_ortho = clipper_core.Coord_orth(*min_coor)
    max_ortho = clipper_core.Coord_orth(*max_coor)
    min_grid = min_ortho.coord_frac(cell).coord_grid(grid)
    max_grid = max_ortho.coord_frac(cell).coord_grid(grid)    
    
def make_unit_cell(model, hklinfo, cell):
    from chimerax.core.geometry import Place, Places
    coord = model.bounds().center().astype(float)
    coord_orth = clipper_core.Coord_orth(*coord)
    coord_frac = coord_orth.coord_frac(cell)
    sg = hklinfo.spacegroup()
    unit_cell_frac_symops = sg.unit_cell_RTops(coord_frac)
    uc_places = []
    for i in range(sg.num_symops()):
        this_op = sg.unit_cell_RTop(unit_cell_frac_symops, i)
        this_op_orth = this_op.rtop_orth(cell)
        uc_places.append(Place(matrix=this_op_orth.matrix()[0:3,:]))
    ucp = Places(uc_places)
    model.positions = ucp
    return ucp
