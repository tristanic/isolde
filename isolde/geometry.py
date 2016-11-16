
# Get the dihedral angle between the vectors p0-p1 and p2-p3 about the
# p1-p2 axis. p0, p1, p2 and p3 should be (x,y,z) coordinates as numpy
# arrays.

    
def get_dihedral(p0, p1, p2, p3):
    import numpy as np
    
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    
    b1 /= np.linalg.norm(b1)
    
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1
    
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.arctan2(y, x)

#def get_dihedral(p0, p1, p2, p3):
    #import numpy as np
    #import warnings
    #warnings.filterwarnings('error')
    #p = np.ones([4,3])
    #for i, pn in enumerate((p0,p1,p2,p3)):
        #p[i] = pn
    #b = p[:-1] - p[1:]
    #b[0] *= -1
    #v = np.array( [ v - (v.dot(b[1])/b[1].dot(b[1])) * b[1] for v in [b[0], b[2]] ] )
    ## Normalize vectors
    #try:
        #v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1,1)
    #except RuntimeWarning:
        #print(p0)
        #print(p1)
        #print(p2)
        #print(p3)
        #raise
    #b1 = b[1] / np.linalg.norm(b[1])
    #x = np.dot(v[0], v[1])
    #m = np.cross(v[0], b1)
    #y = np.dot(m, v[1])
    #return np.arctan2( y, x )
    



#import os
#import ctypes
#_geometry = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)),'./_geometry.so'))

#COORTYPE = ctypes.POINTER(ctypes.c_double * 3)
  
#_get_dihedral = _geometry.get_dihedral
#_get_dihedral.argtypes = [COORTYPE, COORTYPE, COORTYPE, COORTYPE]
#_get_dihedral.restype = ctypes.c_double
   
#def get_dihedral(p0, p1, p2, p3):
    #return _get_dihedral(p0.ctypes.data_as(COORTYPE), 
                    #p1.ctypes.data_as(COORTYPE),
                    #p2.ctypes.data_as(COORTYPE),
                    #p3.ctypes.data_as(COORTYPE))
        
    
