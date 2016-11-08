
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
   
