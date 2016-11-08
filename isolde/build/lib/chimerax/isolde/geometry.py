
# Get the dihedral angle between the vectors p1-p2 and p3-p4 about the
# p2-p3 axis. p1, p2, p3 and p4 should be (x,y,z) coordinates as numpy
# arrays.
def get_dihedral(p0, p1, p2, p3):
    import numpy
    from chimerax.core.geometry import vector as v
    p0p1 = v.normalize_vector(p1 - p0)
    p1p2 = v.normalize_vector(p2 - p1)
    p2p3 = v.normalize_vector(p3 - p2)
    
    norm_1 = v.cross_product(p0p1, p1p2)
    norm_2 = v.cross_product(p1p2, p2p3)
    
    x = v.inner_product(norm_1, norm_2)
    m1 = v.cross_product(norm_1, p1p2)
    y = v.inner_product(m1, norm_2)
    
    import math
    
    return math.atan2(y,x)
    
def get_dihedral_v2(p0, p1, p2, p3):
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

def compare_dihedral_methods(p0, p1, p2, p3, n = 1000):
    from time import time
    start_time = time()
    for i in range(n):
        d = get_dihedral(p0, p1, p2, p3)
    end_time = time()
    print('original: ' +str(d) +' ' + str(end_time - start_time))
    
    start_time = time()
    for i in range(n):
        get_dihedral_v2(p0, p1, p2, p3)
    end_time = time()
    print('original: ' + str(d) + ' ' + str(end_time - start_time))
    
