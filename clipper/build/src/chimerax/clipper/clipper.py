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
    from time import time
    start_time = time()
    n = len(atom_list)
    elements = atom_list.element_names.tolist()
    coords = atom_list.coords
    occupancies = atom_list.occupancy
    u_iso = atom_list.bfactors
    #u_aniso = sel.aniso_u
    u_aniso = numpy.random.rand(n,3,3).astype(float) # FIXME For testing only
    clipper_atom_list = clipper_core.Atom_list(n)

    clipper_atom_list.set_elements(elements)
    clipper_atom_list.set_coord_orth(coords)
    clipper_atom_list.set_occupancies(occupancies)
    clipper_atom_list.set_u_isos(u_iso)
    # need to extract the relevant values into an n x 6 array
    # print('Intermediate time: {}'.format(time() - start_time))
    ua_vals = numpy.empty([n,6], numpy.float)
    ua_vals[:,0] = u_aniso[:,0,0]
    ua_vals[:,1] = u_aniso[:,1,1]
    ua_vals[:,2] = u_aniso[:,2,2]
    ua_vals[:,3] = u_aniso[:,0,1]
    ua_vals[:,4] = u_aniso[:,0,2]
    ua_vals[:,5] = u_aniso[:,1,2]
    clipper_atom_list.set_u_anisos(ua_vals)
    return clipper_atom_list

def map_import_test(filename):
    myhkl = clipper_core.HKL_info()
    fphidata =  clipper_core.HKL_data_F_phi_double()  
    mtzin = clipper_core.CCP4MTZfile()
    mtzin.open_read(filename)
    mtzin.import_hkl_info(myhkl)
    mtzin.import_hkl_data(fphidata, '/crystal/dataset/[2FOFCWT, PH2FOFCWT]')
    mtzin.close_read()
    mygrid = clipper_core.Grid_sampling(myhkl.spacegroup(), myhkl.cell(), myhkl.resolution())
    mymap = clipper_core.Xmap_double(myhkl.spacegroup(), myhkl.cell(), mygrid)
    mymap.fft_from(fphidata)
    return (myhkl, mymap)
    
def vol_box(hklinfo, xmap, min_coor, max_coor):
    cell = hklifno.cell()
    grid = xmap.grid_sampling()
    min_ortho = clipper_core.Coord_orth(*min_coor)
    max_ortho = clipper_core.Coord_orth(*max_coor)
    min_grid = min_ortho.coord_frac(cell).coord_grid(grid)
    max_grid = max_ortho.coord_frac(cell).coord_grid(grid)    
    
