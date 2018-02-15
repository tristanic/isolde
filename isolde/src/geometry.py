# Copyright 2017 Tristan Croll
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


import numpy
import sys, os, glob
import ctypes

NPY_BOOL = ctypes.c_uint8


# Fully pythonic version.
def get_dihedral(p0, p1, p2, p3):
    '''
     Get the dihedral angle between the vectors p0-p1 and p2-p3 about the
    p1-p2 axis. p0, p1, p2 and p3 should be (x,y,z) coordinates as numpy
    arrays.
    '''
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2

    b1 /= numpy.linalg.norm(b1)

    v = b0 - numpy.dot(b0, b1)*b1
    w = b2 - numpy.dot(b2, b1)*b1

    x = numpy.dot(v, w)
    y = numpy.dot(numpy.cross(b1, v), w)
    return numpy.arctan2(y, x)

# C version. Saves about 16 microseconds per dihedral.

dpath = os.path.dirname(os.path.abspath(__file__))
libfile = glob.glob(os.path.join(dpath, '_geometry.cpython*'))[0]
_geometry = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), libfile))
COORTYPE = ctypes.POINTER(ctypes.c_double * 3)

_get_dihedral = _geometry.get_dihedral
_get_dihedral.argtypes = [COORTYPE, COORTYPE, COORTYPE, COORTYPE]
_get_dihedral.restype = ctypes.c_double

def get_dihedral(p0, p1, p2, p3):
    return _get_dihedral(p0.ctypes.data_as(COORTYPE),
                    p1.ctypes.data_as(COORTYPE),
                    p2.ctypes.data_as(COORTYPE),
                    p3.ctypes.data_as(COORTYPE))

# Get values for an entire set of dihedrals at once. More than an order
# of magnitude faster than calling get_dihedral() in a Python loop.
_get_dihedrals = _geometry.get_dihedrals

def get_dihedrals(coords, n):
    full_len = n*4*3
    INTYPE = ctypes.POINTER(ctypes.c_double*full_len)
    ret = numpy.empty(n, numpy.double)
    RTYPE = ctypes.POINTER(ctypes.c_double)
    _get_dihedrals.argtypes = [INTYPE, ctypes.c_int, RTYPE]
    _get_dihedrals(numpy.reshape(coords, -1).ctypes.data_as(INTYPE), n, ret.ctypes.data_as(RTYPE))
    return ret


_rotations = _geometry.rotations
def rotations(axis, angles):
    '''
    Get the rotation matrices around the given axis for an array of angles.
    '''
    n = len(angles)
    ANGLE_TYPE = ctypes.POINTER(ctypes.c_double*n)
    ret = numpy.empty((n,3,4), numpy.double)
    RTYPE=ctypes.POINTER(ctypes.c_double)
    _rotations.argtypes = [COORTYPE, ANGLE_TYPE, ctypes.c_int, RTYPE]
    _rotations(axis.ctypes.data_as(COORTYPE), angles.ctypes.data_as(ANGLE_TYPE),
               n, ret.ctypes.data_as(RTYPE))
    return ret

_scale_transforms = _geometry.scale_transforms
def scale_transforms(scales, transforms):
    ''' Get a set of transformations scaling coordinates by the values in scales.'''
    n = len(scales)
    tf = convert_and_sanitize_numpy_array(transforms, numpy.double)
    SCALE_TYPE = ctypes.POINTER(ctypes.c_double*n)
    RTYPE = ctypes.POINTER(ctypes.c_double)
    _scale_transforms.argtypes = [SCALE_TYPE, ctypes.c_int, RTYPE]
    _scale_transforms(scales.ctypes.data_as(SCALE_TYPE), n, tf.ctypes.data_as(RTYPE))
    return tf



_multiply_transforms = _geometry.multiply_transforms
def multiply_transforms(tf1, tf2):
    TF_TYPE = ctypes.POINTER(ctypes.c_double*12)
    ret = numpy.empty((3,4), numpy.double)
    _multiply_transforms.argtypes = [TF_TYPE, TF_TYPE, TF_TYPE]
    _multiply_transforms(tf1.ctypes.data_as(TF_TYPE),
                         tf2.ctypes.data_as(TF_TYPE),
                         ret.ctypes.data_as(TF_TYPE))
    return ret

def convert_and_sanitize_numpy_array(array, dtype):
    '''
    Convert a numpy array to the specified data type, and ensure its
    contents are C-contiguous in memory.
    '''
    #~ if array.flags.c_contiguous:
        #~ if array.dtype == dtype:
            #~ return array
        #~ return array.as_type(dtype)
    ret = numpy.empty(array.shape, dtype)
    ret[:] = array
    return ret


_flip_rotate_shift = _geometry.flip_rotate_and_shift
def flip_rotate_and_shift(flip_mask, flip_tf, rotations, shifts):
    n = len(flip_mask)
    if len(rotations) != n or len(shifts) != n:
        raise TypeError('flip_mask, rotations and shifts must all be the '\
            +'same length!')
    #~ ftf = numpy.empty(flip_tf.shape, numpy.double)
    #~ ftf[:] = flip_tf
    from numpy import double
    ftf = convert_and_sanitize_numpy_array(flip_tf, double)
    rot = convert_and_sanitize_numpy_array(rotations, double)
    sh = convert_and_sanitize_numpy_array(shifts, double)
    #~ shape = rotations.shape
    #~ rot = numpy.empty(shape, numpy.double)
    #~ rot[:] = rotations
    #~ sh = numpy.empty(shape, numpy.double)
    #~ sh[:] = shifts
    FLAG_TYPE = ctypes.POINTER(NPY_BOOL*n)
    TF_TYPE = ctypes.POINTER(ctypes.c_double*12)
    TF_ARRAY_TYPE = ctypes.POINTER(ctypes.c_double*12*n)
    ret = numpy.empty(rotations.shape, numpy.double)
    _flip_rotate_shift.argtypes = [ctypes.c_int, FLAG_TYPE, TF_TYPE,
                                   TF_ARRAY_TYPE, TF_ARRAY_TYPE, TF_ARRAY_TYPE]
    _flip_rotate_shift(n, flip_mask.ctypes.data_as(FLAG_TYPE),
                          ftf.ctypes.data_as(TF_TYPE),
                          rot.ctypes.data_as(TF_ARRAY_TYPE),
                          sh.ctypes.data_as(TF_ARRAY_TYPE),
                          ret.ctypes.data_as(TF_ARRAY_TYPE))
    return ret

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

    n = cross(varray[1]-varray[0], varray[3]-varray[0])

    narray[:] = n
    tarray = array([[0,1,4],[1,2,4],[2,3,4]],int32)
    return varray, narray, tarray


def dihedral_fill_planes(dihedrals, target_drawing):
    dw = target_drawing
    varray = numpy.empty([5*len(dihedrals), 3], numpy.float32)
    narray = numpy.empty([5*len(dihedrals), 3], numpy.float32)
    tarray = numpy.empty([3*len(dihedrals), 3], numpy.int32)
    ntri = 0
    for i, d in enumerate(dihedrals):
        thisv, thisn, thist = dihedral_fill_plane(*d.coords)
        start, end = i*5, i*5+5
        varray[start:end] = thisv
        narray[start:end] = thisn
        thist += i*5
        tarray[ntri:ntri+3] = thist
        ntri += 3
    dw.vertices, dw.normals, dw.triangles = varray, narray, tarray

_dihedral_fill_planes=_geometry.dihedral_fill_planes
def dihedral_fill_planes(dihedrals, target_drawing):
    dw = target_drawing
    coords = dihedrals.coords
    n = len(dihedrals)
    varray = numpy.empty([5*n, 3], numpy.float32)
    narray = numpy.empty([5*n, 3], numpy.float32)
    tarray = numpy.empty([3*n, 3], numpy.int32)
    COORTYPE = ctypes.POINTER(ctypes.c_double*coords.size)
    V_TYPE = ctypes.POINTER(ctypes.c_float*varray.size)
    N_TYPE = ctypes.POINTER(ctypes.c_float*narray.size)
    T_TYPE = ctypes.POINTER(ctypes.c_int32*tarray.size)
    _dihedral_fill_planes.argtypes = [ctypes.c_int, COORTYPE, V_TYPE, N_TYPE, T_TYPE]
    _dihedral_fill_planes(n, coords.ctypes.data_as(COORTYPE),
                             varray.ctypes.data_as(V_TYPE),
                             narray.ctypes.data_as(N_TYPE),
                             tarray.ctypes.data_as(T_TYPE))
    dw.vertices, dw.normals, dw.triangles = varray, narray, tarray

_dihedral_fill_and_color_planes = _geometry.dihedral_fill_and_color_planes
def dihedral_fill_and_color_planes(dihedrals, target_drawing,
            twisted_mask, cis_pro_mask, cis_nonpro_color, twisted_color,
            cis_pro_color):
    dw = target_drawing
    coords = dihedrals.coords
    n = len(dihedrals)
    varray = numpy.empty([5*n, 3], numpy.float32)
    narray = numpy.empty([5*n, 3], numpy.float32)
    tarray = numpy.empty([3*n, 3], numpy.int32)
    carray = numpy.empty([5*n, 4], numpy.uint8)
    COORTYPE = ctypes.POINTER(ctypes.c_double*coords.size)
    V_TYPE = ctypes.POINTER(ctypes.c_float*varray.size)
    N_TYPE = ctypes.POINTER(ctypes.c_float*narray.size)
    T_TYPE = ctypes.POINTER(ctypes.c_int32*tarray.size)
    C_TYPE = ctypes.POINTER(ctypes.c_uint8*carray.size)
    COLOR_TYPE = ctypes.POINTER(ctypes.c_uint8*4)
    MASK_TYPE = ctypes.POINTER(NPY_BOOL*n)
    _dihedral_fill_planes.argtypes = [ctypes.c_int, COORTYPE, MASK_TYPE,
        MASK_TYPE, COLOR_TYPE, COLOR_TYPE, COLOR_TYPE,
        V_TYPE, N_TYPE, T_TYPE, C_TYPE]
    _dihedral_fill_and_color_planes(n, coords.ctypes.data_as(COORTYPE),
                             twisted_mask.ctypes.data_as(MASK_TYPE),
                             cis_pro_mask.ctypes.data_as(MASK_TYPE),
                             cis_nonpro_color.ctypes.data_as(COLOR_TYPE),
                             twisted_color.ctypes.data_as(COLOR_TYPE),
                             cis_pro_color.ctypes.data_as(COLOR_TYPE),
                             varray.ctypes.data_as(V_TYPE),
                             narray.ctypes.data_as(N_TYPE),
                             tarray.ctypes.data_as(T_TYPE),
                             carray.ctypes.data_as(C_TYPE))
    dw.vertices, dw.normals, dw.triangles, dw.vertex_colors = varray, narray, tarray, carray

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

def exclamation_mark(radius = 0.1, height = 2, nc=8, color = [255,0,0,255]):
    '''
    An exclamation mark for annotating rotamer outliers
    '''
    from chimerax.core.surface.shapes import cone_geometry, sphere_geometry2
    from chimerax.core.geometry import translation, scale
    import numpy
    stem = cone_geometry(radius=radius, height=height, nc=nc, caps=True)
    spheres = list(sphere_geometry2(nc*4))
    spheres[0] = scale(radius*0.7).moved(spheres[0])

    v, n, t = stem


    vbottom = translation((0,0, height/2+radius*1.5)).moved(spheres[0])
    t = numpy.concatenate((t, spheres[2]+len(v)))
    v = numpy.concatenate((v, vbottom))
    n = numpy.concatenate((n, spheres[1]))

    return v, n, t


def spiral(major_radius=0.25, minor_radius=0.1, height=1, turns=1.0,
           turn_segments=10, circle_segments=5):
    '''
    Draw a 3D spiral.
    '''
    from math import pi
    from chimerax.core.surface import tube
    spline_points = numpy.array(range(turn_segments), float)/(turn_segments-1)
    spline_xyz = numpy.empty((turn_segments, 3), float)
    spline_xyz[:,0] = numpy.cos(spline_points * turns * 2*pi)*major_radius
    spline_xyz[:,1] = numpy.sin(spline_points * turns * 2*pi)*major_radius
    spline_xyz[:,2] = spline_points * height
    return tube.tube_spline(spline_xyz, minor_radius, segment_subdivisions = 2, circle_subdivisions = circle_segments)


def simple_arrow(radius = 0.1, height = 1, nc = 20, color = [255, 0, 0, 255], caps = True,
                    head_length_fraction = 0.33, head_width_ratio = 1.5,
                    points_out = True):
    from chimerax.core.models import Drawing
    d = Drawing(name='Arrow')
    d.color = color
    d.vertices, d.normals, d.triangles = simple_arrow_geometry(
        radius, height, nc, caps, head_length_fraction,
        head_width_ratio, points_out)
    return d

def simple_arrow_geometry(radius = 0.1, height = 1, nc = 20, caps = True,
                    head_length_fraction = 0.33, head_width_ratio = 1.5,
                    points_out = True):
    '''
    Define a simple 3D arrow made from two cones joined base-to-base.
    If points_out is true the point of the longer, narrower cone will
    be at the origin.
    '''
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
    #~ head.vertices = hver
    #~ head.normals = hn
    #~ head.triangles = ht
    #~ shaft.vertices = sver
    #~ shaft.normals = sn
    #~ shaft.triangles = st
    #~ d.add_drawing(head)
    #~ d.add_drawing(shaft)
    v = numpy.concatenate((hver, sver))
    n = numpy.concatenate((hn, sn))
    t = numpy.concatenate((ht, st+len(hver)))

    return v,n,t

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

def pin_geometry(handle_radius, pin_radius, total_height):
    '''
    Simple 3D representation of a drawing pin.
    Args:
        handle_radius:
            The radius of the "handle" in Angstroms
        pin_radius:
            The radius of the "pin" in Angstroms
        height:
            The overall height of the drawing
    '''
    import numpy
    from chimerax.core.surface.shapes import cylinder_geometry, cone_geometry
    from chimerax.core.geometry import translation
    pin_height = total_height*5/12
    handle_height = total_height*7/12
    tb_height = total_height/4

    pin = list(cone_geometry(radius = pin_radius, height = pin_height, points_up = False))
    handle_bottom = list(cone_geometry(radius = handle_radius, height = tb_height, nc = 8))
    handle_middle = list(cylinder_geometry(radius = handle_radius/2, height = handle_height, nc = 8, caps = False))
    handle_top = list(cone_geometry(radius = handle_radius, height = tb_height, nc = 8, points_up = False))


    pint = translation((0,0,pin_height/2))
    pin[0] = pint.moved(pin[0])
    hbt = translation((0,0, pin_height + tb_height/2.05))
    handle_bottom[0] = hbt.moved(handle_bottom[0])
    hmt = translation((0,0, handle_height/2 + pin_height))
    handle_middle[0] = hmt.moved(handle_middle[0])
    htt = translation((0,0,total_height - tb_height/2.05))
    handle_top[0] = htt.moved(handle_top[0])

    vertices = numpy.concatenate((pin[0], handle_bottom[0], handle_middle[0], handle_top[0]))
    normals = numpy.concatenate((pin[1], handle_bottom[1], handle_middle[1], handle_top[1]))

    ntri = len(pin[0])
    triangles = pin[2]
    for d in (handle_bottom, handle_middle, handle_top):
        triangles = numpy.concatenate((triangles, d[2]+ntri))
        ntri += len(d[0])

    return vertices, normals, triangles

def dumbbell_geometry(major_radius=1, minor_radius=0.1, thickness=0.1, height=1, nz=2, nc1 = 10, nc2=10):
    from chimerax.core.surface.shapes import cylinder_geometry
    dva, dna, dta = cylinder_geometry(major_radius, thickness, nz, nc1, True)
    hva, hna, hta = cylinder_geometry(minor_radius, height, nz, nc2, False)
    nv = len(dva)
    vs = []
    ns = []
    ts = []
    offset = numpy.array([0,0,height/2])
    vs.append(dva-offset)
    ns.append(dna)
    ts.append(dta)
    vs.append(hva)
    ns.append(hna)
    ts.append(hta + nv)
    nv += len(hva)
    vs.append(dva+offset)
    ns.append(dna)
    ts.append(dta+nv)
    vd, nd, td = numpy.concatenate(vs), numpy.concatenate(ns), numpy.concatenate(ts)
    return vd, nd, td

def arc_points(n, radius, starting_angle, final_angle):
    final_angle = final_angle*n/(n-1)
    from numpy import arange, float32, empty, sin, cos
    from math import pi
    a = arange(n) * ((final_angle-starting_angle)/n) + starting_angle
    c = empty((n,3), float32)
    c[:,0] = radius*cos(a)
    c[:,1] = radius*sin(a)
    c[:,2] = 0
    return c


def split_torus_geometry(major_radius, minor_radius, circle_segments, ring_segments, starting_angle, final_angle):
    '''
    Define a torus (doughnut), where major_radius is the distance from
    the centre of the ring to the axis of the solid portion, and
    minor_radius defines the thickness of the solid_portion.
    '''
    from math import pi
    from chimerax.core.surface import tube
    from chimerax.core.geometry import rotation
    path = arc_points(ring_segments, major_radius, starting_angle, final_angle)
    return tube.tube_spline(path, minor_radius, segment_subdivisions = 2, circle_subdivisions = circle_segments)

def ring_arrow(major_radius, minor_radius, circle_segments, ring_segments, head_length, head_radius):
    import numpy
    #Find the starting angle from the length of the arrow head
    from math import asin, pi, degrees
    starting_angle = 2*asin(head_length/(2*major_radius))
    vr, nr, tr = split_torus_geometry(major_radius, minor_radius, circle_segments, ring_segments, starting_angle*0.8, 3*pi/2)
    vh, nh, th = cone_geometry(head_radius, head_length, circle_segments)
    from chimerax.core.geometry import rotation
    # Rotate the cone
    move_dir = numpy.array(((major_radius, 0, -head_length),), numpy.double)
    r1 = rotation((1,0,0), -90)
    r1.move(vh)
    r1.move(move_dir)
    r2 = rotation((0,0,1), degrees(starting_angle))
    r2.move(vh)
    r2.move(move_dir)
    # ... and move it into position
    vh += move_dir

    from numpy import concatenate
    v = concatenate((vr, vh))
    n = concatenate((nr, nh))
    t = concatenate((tr, th+len(vr)))

    return v, n, t

def ring_arrow_with_post(major_radius, minor_radius, circle_segments,
                         ring_segments, head_length, head_radius,
                         post_radius, post_height):
    v, n, t = ring_arrow(major_radius, minor_radius, circle_segments,
                         ring_segments, head_length, head_radius)
    pv, pn, pt = post_geometry(post_radius, post_height, caps=True)
    from numpy import concatenate
    rv = concatenate((v, pv))
    rn = concatenate((n, pn))
    rt = concatenate((t, pt+len(v)))
    return rv, rn, rt

def test_ra(session):
    from chimerax.core.models import Model, Drawing
    m = Model('test', session)
    d = Drawing('ring')
    d.vertices, d.normals, d.triangles = ring_arrow_with_post(0.5, 0.05, 4, 6, 0.3, 0.1, 0.05, 1)
    m.add_drawing(d)
    session.models.add([m])

def post_geometry(radius, height, caps=False):
    '''
    Returns a simple cylinder, rotated and translated so its base is on
    the origin and it points along (1,0,0)
    '''
    from chimerax.core.surface.shapes import cylinder_geometry
    from chimerax.core.geometry import rotation, translation
    v, n, t = cylinder_geometry(radius=radius, height=height, caps=caps, nc=6)
    tr = translation([0,0,height/2])
    tr.move(v)
    r = rotation([0,1,0], 90)
    r.move(v)
    r.apply_without_translation(n)
    return v,n,t

def bond_cylinder_placements(bonds):
    '''From chimerax.core.structure._bond_cylinder_placements.'''

    n = len(bonds)
    from numpy import empty, float32
    p = empty((n,4,4), float32)

    radii = numpy.ones(len(bonds))
    from chimerax.core.geometry import cylinder_rotations, Places
    axyz0, axyz1 = [a.coords for a in bonds.atoms]
    cylinder_rotations(axyz0, axyz1, radii, p)

    p[:,3,:3] = 0.5*(axyz0 + axyz1)

    pl = Places(opengl_array = p)
    return pl
