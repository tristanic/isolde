

# Fully pythonic version.    
def get_dihedral(p0, p1, p2, p3):
    '''
     Get the dihedral angle between the vectors p0-p1 and p2-p3 about the
    p1-p2 axis. p0, p1, p2 and p3 should be (x,y,z) coordinates as numpy
    arrays.
    '''
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

# C version. Saves about 16 microseconds per dihedral.
#import os
#import ctypes
#_geometry = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)),'./lib_geometry.so'))

#COORTYPE = ctypes.POINTER(ctypes.c_double * 3)
  
#_get_dihedral = _geometry.get_dihedral
#_get_dihedral.argtypes = [COORTYPE, COORTYPE, COORTYPE, COORTYPE]
#_get_dihedral.restype = ctypes.c_double
   
#def get_dihedral(p0, p1, p2, p3):
    #return _get_dihedral(p0.ctypes.data_as(COORTYPE), 
                    #p1.ctypes.data_as(COORTYPE),
                    #p2.ctypes.data_as(COORTYPE),
                    #p3.ctypes.data_as(COORTYPE))
        
def dihedral_fill_plane(p0, p1, p2, p3):
    '''
    Fill in the "cup" in a dihedral with a pseudo-planar surface
    '''
    from numpy import empty, float32, int32, array, cross
    varray = empty((5,3), float32)
    narray = empty((5,3), float32)
    
    for i, p in enumerate((p0,p1,p2,p3)):
        varray[i] = p
    
    varray[4] = (varray[0]+varray[3])/2
    
    # This surface should always be almost planar, so we'll assign it a
    # single normal
    
    n = cross(varray[3]-varray[0], varray[1]-varray[0])
    
    narray[:] = n
    tarray = array([[0,1,4],[1,2,4],[2,3,4]],int32)
    return varray, varray, tarray
    



def cone_geometry(radius = 1, height = 1, nc = 10, caps = True, flipped = False):
    '''
    Return vertex, normal vector and triangle arrays for cone geometry
    with specified radius and height. If flipped is true, the base of the
    cone will be at the origin, otherwise the point will be at the origin.
    '''
    from numpy import ones, empty, float32, arange, cos, sin, int32, pi
    vc = nc * 2
    tc = nc
    if caps:
        vc += (nc + 1)
        tc += nc
    varray = empty((vc, 3), float32)
    narray = empty((vc, 3), float32)
    tarray = empty((tc, 3), int32)

    # Compute a circle (which may be used twice if caps is true)
    angles = (2 * pi / nc) * arange(nc)
    import sys
    circle = empty((nc, 2), float32)
    circle[:,0] = cos(angles) * radius
    circle[:,1] = sin(angles) * radius

    # Create cone faces (first nc*2 vertices)
    # The normals are wrong, but let's see how they look
    nc2 = nc * 2
    if not flipped:
        varray[:nc] = (0, 0, 0)      # point of cone (multiple normals)
    else:
        varray[:nc] = (0, 0, height)
    narray[:nc,:2] = circle
    if not flipped:
        narray[:nc,2] = 0
    else:
        narray[:nc,2] = height            
    varray[nc:nc2,:2] = circle      # base of cone
    if not flipped:
        varray[nc:nc2,2] = height
    else:
        varray[nc:nc2,2] = 0
    narray[nc:nc2,:2] = circle      # wrong, but close (enough?)
    if not flipped:
        narray[nc:nc2,2] = height
    else:
        narray[nc:nc2,2] = 0
    tarray[:nc,0] = arange(nc)
    tarray[:nc,1] = (arange(nc) + 1) % nc + nc
    tarray[:nc,2] = arange(nc) + nc

    # Create cone base (last nc+1 vertices)
    if caps:
        if not flipped:
            varray[nc2] = (0, 0, height)
        else:
            varray[nc2] = (0, 0, 0)
        varray[nc2+1:,:2] = circle
        if not flipped:
            varray[nc2+1:,2] = height
        else:
            varray[nc2+1:,2] = 0
        narray[nc2:] = (0, 0, 1)
        tarray[nc:,0] = nc2
        tarray[nc:,1] = (arange(nc) + 1) % nc + nc2 + 1
        tarray[nc:,2] = arange(nc) + nc2 + 1

    return varray, narray, tarray
    
def simple_arrow(radius = 0.1, height = 1, nc = 20, color = [255, 0, 0, 255], caps = True, 
                    head_length_fraction = 0.33, head_width_ratio = 1.5, 
                    points_out = True):
    '''
    Define a simple 3D arrow made from two cones joined base-to-base.
    If points_out is true the point of the longer, narrower cone will
    be at the origin.
    '''
    from chimerax.core.models import Drawing
    d = Drawing(name='Arrow')
    head = Drawing(name = 'arrow_head')
    shaft = Drawing(name = 'arrow_shaft')
    head.color = color
    shaft.color = color
    head_base_width = radius
    shaft_base_width = radius / head_width_ratio
    head_length = height * head_length_fraction
    shaft_length = height * (1 - head_length_fraction)
    
    nc2 = nc * 2
    if points_out:
        hver, hn, ht = cone_geometry(head_base_width, head_length, nc, flipped = True)
        sver, sn, st = cone_geometry(shaft_base_width, shaft_length, nc, flipped = False)
        # shift the head outwards
        hver[:nc] = (0,0,height)
        hn[:nc,:2] = height
        hver[nc:nc2,2] = shaft_length
        hn[nc:nc2,2] = shaft_length
        hver[nc2] = (0, 0, shaft_length)
        hver[nc2+1:,2] = shaft_length
    else:
        hver, hn, ht = cone_geometry(head_base_width, head_length, nc, flipped = False)
        sver, sn, st = cone_geometry(shaft_base_width, shaft_length, nc, flipped = True)
        # shift the head outwards
        sver[:nc] = (0,0,height)
        sn[:nc,:2] = height
        sver[nc:nc2,2] = head_length
        sn[nc:nc2,2] = head_length
        sver[nc2] = (0, 0, head_length)
        sver[nc2+1:,2] = head_length
        
    import numpy
    head.vertices = hver
    head.normals = hn
    head.triangles = ht
    shaft.vertices = sver
    shaft.normals = sn
    shaft.triangles = st
    d.add_drawing(head)
    d.add_drawing(shaft)
    #d.vertices = numpy.concatenate((hver, sver))
    #d.normals = numpy.concatenate((hn, sn))
    #d.triangles = numpy.concatenate((ht, st))
    
    return d
    
def arrow_between_points(arrow, xyz0, xyz1):
    '''
    Takes an arrow drawn by simple_arrow and rotates, translates and scales
    it to point from xyz0 to xyz1 (or vice versa, if it's an inward-pointing
    arrow). In either case the arrow will have one end on xyz0. The other end
    will be on xyz1 if the height of the original arrow was 1 unit.
    '''
    from chimerax.core.geometry.place import translation, vector_rotation, scale, Place, product
    import numpy
    
    # Original arrow points along the z axis
    u = numpy.array((0,0,1),dtype='float32')
    # Get vector from xyz0 to xyz1
    v = xyz1-xyz0
    # Get vector length
    l = numpy.linalg.norm(v)
    rot = vector_rotation(u,v/l)
    trans = translation(xyz0)
    sc = scale(l)
    arrow.position = product((trans,rot,sc))
    
def arrow_along_force_vector(arrow, xyz0, force, scale = 1.0):
    '''
    Takes an arrow drawn by simple_arrow and rotates, translates and scales
    it to point from xyz0 along a force vector, with length scaled according
    to the magnitude of the force.
    '''
    from chimerax.core.geometry.place import translation, vector_rotation, scale, Place, product
    import numpy
    
    # Original arrow points along the z axis
    u = numpy.array((0,0,1),dtype='float32')
    # Get vector length
    l = numpy.linalg.norm(force)
    rot = vector_rotation(u,force/l)
    trans = translation(xyz0)
    sc = scale(l)
    arrow.position = product((trans,rot,sc))
    
    
        
