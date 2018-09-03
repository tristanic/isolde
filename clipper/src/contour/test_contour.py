def test_contour(session, level=0.2, copy=False):
    import numpy
    from chimerax.clipper.symmetry import get_symmetry_handler
    sh = get_symmetry_handler(session.isolde.selected_model)
    v = sh.xmapset[1]
    from chimerax.clipper.contour_thread import Contour_Thread_Mgr
    cm = Contour_Thread_Mgr()
    matrix = v.matrix()
    m = None
    if matrix.dtype != numpy.float32:
        m = numpy.empty(matrix.shape, numpy.float32)
        m[:] = matrix
    else:
        m = matrix
    cm.start_compute(m, level, False, True)
    va, ta, na = cm.get_result()
    #return va, ta, na
    va, na, ta, hidden_edges = v.surfaces[0]._adjust_surface_geometry(va, na, ta, v.rendering_options, 0.233)
    v.surfaces[0]._set_surface(va, na, ta, hidden_edges)
    return va, ta, na
