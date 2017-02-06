from . import clipper
import numpy


# Override logging of Clipper messages to send them to the ChimeraX log
# rather than just printing them to the console.

_clipper_messages = clipper._clipper_messages

def _log_clipper(func):
    def func_wrapper(*args, **kwargs):
        _clipper_messages.clear()
        func(*args, **kwargs)
        message_string = _clipper_messages.read_and_clear()
        if message_string:
            session.logger.info("CLIPPER WARNING:")
            session.logger.info(message_string)
    return func_wrapper

clipper.log_clipper = _log_clipper

'''
class Clipper_Project:
    def __init__(self, project_name):
        from .data_tree import DataTree
        self.data = Datatree()
        self.data['Head'].set_metadata(('Name', 'project_name'))
'''        
    
from .data_tree import DataTree
class Clipper_MTZ(DataTree):
    '''
    A container class to capture all the information loaded in from a 
    .mtz or .cif structure factor file.
    structure, using the clipper-python libraries: 
        - The reflection data
        - The cell and symmetry parameters
        - The map(s)
        - The atomic model describing the asymmetric unit
    
    Contains methods for:
        - Loading in reflection and reciprocal-space map data from 
          structure factor (.mtz or .cif) files
        - Displaying density surrounding the current center of rotation
        - Quickly finding and displaying all local symmetry equivalents
        ...
    '''
    
    
     
    
    
    
    def __init__(self):
        DataTree.__init__(self)
        
    def load_hkl_data(session, filename):
        import os
        if os.path.isfile(filename):
            hklfile = os.path.abspath(filename)
        else:
            self.session.logger.info('Invalid filename!')
            return
        
        extension = os.path.splitext(filename)[1].lower()
        if extension not in ('.cif', '.mtz'):
            self.session.logger.info('Reflection file must be in either .cif or .mtz format!')
            return
        
        from datetime import datetime
        self.set_metadata(('Source file', hklfile))
        self.set_metadata(('Accessed',str(datetime.today())))
        
        if extension == '.mtz':
            mtzin = clipper.CCP4MTZfile()
            hkl = self['hkl'] = clipper.HKL_info()
            mtzin.open_read(filename)
            mtzin.import_hkl_info(hkl)
            # Get all the column names and types
            column_labels = mtzin.column_labels
            # Sort the columns into groups, and organise into a temporary tree 
            temp_tree = DataTree()
            i = 0
            for l in column_labels:
                thisname, thistype = l.__str__().split(' ')
                crystal, dataset, name = thisname.split('/')[1:]
                # The h, k, l indices are already captured in self.hklinfo
                if thistype != 'H':
                    temp_tree[crystal][dataset][thistype][name] = i
                    i += 1
            
            for crystal in temp_tree.keys():
                for dataset in temp_tree[crystal].keys():
                    # Find I/sigI pairs and pull them in as
                    # Clipper HKL_data_I_sigI objects
                    for iname in temp_tree[crystal][dataset]['J'].keys():
                        for sname in temp_tree[crystal][dataset]['Q'].keys():
                            if iname in sname:
                                this_isigi \
                                    = self[crystal][dataset]['Intensities'][iname + ', ' + sname] \
                                    = clipper.HKL_data_I_sigI()
                                mtzin.import_hkl_data(this_isigi,
                                '/' + crystal + '/' + dataset + 
                                '/[' + iname + ',' + sname + ']')
                    
                    
                    
                    # Find amplitude/phase pairs and pull them in as 
                    # Clipper HKL_data_F_phi objects
                    for fname in temp_tree[crystal][dataset]['F'].keys():
                        for pname in temp_tree[crystal][dataset]['P'].keys():
                            if fname in pname:
                                thisfphi \
                                    = self[crystal][dataset]['F_Phi'][fname + ', ' + pname] \
                                    = clipper.HKL_data_F_phi()
                                mtzin.import_hkl_data(thisfphi, 
                                '/' + crystal + '/' + dataset + 
                                '/[' + fname + ',' + pname + ']')
                    
            
            
            
            
            
            
            
            
                
            
                
        
        
    
    def import_Xmap_from_mtz_test(session, filename):
        myhkl = HKL_info(session)
        fphidata =  HKL_data_F_phi()  
        mtzin = CCP4MTZfile()
        mtzin.open_read(filename)
        mtzin.import_hkl_info(myhkl)
        mtzin.import_hkl_data(fphidata, '/crystal/dataset/[2FOFCWT, PH2FOFCWT]')
        mtzin.close_read()
        name = '2FOFCWT'
        mygrid = Grid_sampling(myhkl.spacegroup(), myhkl.cell(), myhkl.resolution())
        mymap = Xmap(session, name, myhkl.spacegroup(), myhkl.cell(), mygrid)
        mymap.fft_from(fphidata)
        return (fphidata, myhkl, mymap)
    
    
