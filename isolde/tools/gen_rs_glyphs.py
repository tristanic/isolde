# @Author: Tristan Croll
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

'''
One-off OFFLINE generator for the extruded "R"/"S" letter meshes used by
ISOLDE's chiral R/S markup (``src/validation/_rs_glyphs.py``).

This is NOT shipped or run inside ChimeraX. ChimeraX's bundled Python has no
font->outline facility (no freetype/fontTools; PIL/Qt only rasterise), so the two
letter meshes are precomputed here once and committed as hardcoded numpy arrays
-- giving the bundle real depth-tested 3D letters with zero runtime dependencies.

Run in a throwaway venv with ``fonttools`` + ``numpy`` installed:

    python -m venv /tmp/glyphgen
    /tmp/glyphgen/Scripts/python -m pip install fonttools numpy
    /tmp/glyphgen/Scripts/python tools/gen_rs_glyphs.py

It extracts the 'R' and 'S' glyph outlines from a TrueType font, flattens the
quadratic beziers, triangulates the filled face (hole-bridged ear clipping -- the
'R' counter is a hole), extrudes to a prism (front + back faces + side walls),
normalises each letter to unit height centred at the origin in the XY plane, and
writes the data module. Each letter is validated by asserting the triangulated
front-face area equals the signed polygon area (outer - holes).
'''

import os
import math
import numpy

FONT_PATH = r'C:\Windows\Fonts\arialbd.ttf'   # Arial Bold: clean, legible counters
BEZIER_STEPS = 3         # samples per quadratic segment (small markup glyphs)
DEPTH = 0.18             # extrusion depth in unit-height units
OUT_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    'src', 'validation', '_rs_glyphs.py')


# ---------------------------------------------------------------------------
# Outline extraction (fontTools)
# ---------------------------------------------------------------------------

def glyph_contours(letter):
    '''List of contours (each a list of (x,y)) for ``letter`` in FONT_PATH,
    with quadratic beziers flattened to polylines.'''
    from fontTools.ttLib import TTFont
    from fontTools.pens.recordingPen import RecordingPen
    font = TTFont(FONT_PATH)
    gname = font.getBestCmap()[ord(letter)]
    pen = RecordingPen()
    font.getGlyphSet()[gname].draw(pen)
    contours = []
    cur = []
    last = None
    for op, args in pen.value:
        if op == 'moveTo':
            cur = [args[0]]
            last = args[0]
        elif op == 'lineTo':
            cur.append(args[0])
            last = args[0]
        elif op == 'qCurveTo':
            # TrueType quadratic run: all but the last point are off-curve control
            # points; consecutive off-curve points imply an on-curve midpoint.
            pts = list(args)
            on_end = pts[-1]
            offs = pts[:-1]
            start = last
            # Expand implied on-curve points between consecutive control points.
            segs = []
            for i, c in enumerate(offs):
                if i < len(offs) - 1:
                    nxt = offs[i + 1]
                    mid = ((c[0] + nxt[0]) / 2.0, (c[1] + nxt[1]) / 2.0)
                    segs.append((c, mid))
                else:
                    segs.append((c, on_end))
            for ctrl, end in segs:
                for s in range(1, BEZIER_STEPS + 1):
                    t = s / BEZIER_STEPS
                    mt = 1.0 - t
                    x = mt * mt * start[0] + 2 * mt * t * ctrl[0] + t * t * end[0]
                    y = mt * mt * start[1] + 2 * mt * t * ctrl[1] + t * t * end[1]
                    cur.append((x, y))
                start = end
            last = on_end
        elif op == 'closePath':
            if cur:
                # Drop a duplicated closing point if present.
                if len(cur) > 1 and _close(cur[0], cur[-1]):
                    cur = cur[:-1]
                contours.append(cur)
            cur = []
    return contours


def _close(a, b, eps=1e-6):
    return abs(a[0] - b[0]) < eps and abs(a[1] - b[1]) < eps


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def signed_area(poly):
    s = 0.0
    n = len(poly)
    for i in range(n):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % n]
        s += x0 * y1 - x1 * y0
    return 0.5 * s


def point_in_poly(p, poly):
    x, y = p
    inside = False
    n = len(poly)
    j = n - 1
    for i in range(n):
        xi, yi = poly[i]
        xj, yj = poly[j]
        if ((yi > y) != (yj > y)) and \
                (x < (xj - xi) * (y - yi) / (yj - yi + 1e-30) + xi):
            inside = not inside
        j = i
    return inside


def orient(poly, ccw=True):
    a = signed_area(poly)
    if (a < 0 and ccw) or (a > 0 and not ccw):
        return poly[::-1]
    return poly


# ---------------------------------------------------------------------------
# Hole bridging (Eberly) + ear clipping
# ---------------------------------------------------------------------------

def bridge_holes(outer, holes):
    '''Splice each hole into ``outer`` via a mutually-visible bridge, returning a
    single simple polygon. ``outer`` CCW, ``holes`` CW.'''
    poly = list(outer)
    # Process holes by descending max-x so inner bridges don't cross later ones.
    for hole in sorted(holes, key=lambda h: max(p[0] for p in h), reverse=True):
        poly = _bridge_one(poly, hole)
    return poly


def _bridge_one(outer, hole):
    # Eberly "Triangulation by Ear Clipping", hole-visibility step.
    mi = max(range(len(hole)), key=lambda i: hole[i][0])
    M = hole[mi]
    # Cast +x ray from M; find nearest intersection with an outer edge.
    best_x = math.inf
    best_edge = None
    for i in range(len(outer)):
        a = outer[i]
        b = outer[(i + 1) % len(outer)]
        if (a[1] > M[1]) == (b[1] > M[1]):
            continue
        t = (M[1] - a[1]) / (b[1] - a[1])
        ix = a[0] + t * (b[0] - a[0])
        if ix >= M[0] and ix < best_x:
            best_x = ix
            best_edge = (i, a, b)
    if best_edge is None:
        raise RuntimeError('hole bridge: no visible outer edge')
    i, a, b = best_edge
    I = (best_x, M[1])
    # Candidate P: the edge endpoint with larger x (the visible reflex apex).
    P_idx = i if a[0] > b[0] else (i + 1) % len(outer)
    P = outer[P_idx]
    # Refine: any outer vertex inside triangle (M, I, P) and reflex becomes the
    # better bridge target (closest by angle to the +x ray).
    tri = [M, I, P]
    best_ang = math.inf
    for j, v in enumerate(outer):
        if v in (P,):
            continue
        if point_in_poly(v, tri):
            ang = abs(math.atan2(v[1] - M[1], v[0] - M[0]))
            if ang < best_ang:
                best_ang = ang
                P_idx = j
                P = v
    # Splice: outer[..P] + hole rotated to start at M + M..P bridge + outer[P..].
    hole_seq = hole[mi:] + hole[:mi] + [hole[mi]]
    return outer[:P_idx + 1] + hole_seq + [P] + outer[P_idx + 1:]


def ear_clip(poly):
    '''Triangulate a simple CCW polygon; returns list of (i,j,k) into ``poly``.'''
    n = len(poly)
    idx = list(range(n))
    tris = []
    guard = 0
    while len(idx) > 3 and guard < 10000:
        guard += 1
        ear = False
        m = len(idx)
        for k in range(m):
            i0, i1, i2 = idx[(k - 1) % m], idx[k], idx[(k + 1) % m]
            a, b, c = poly[i0], poly[i1], poly[i2]
            # Convex vertex? (CCW polygon -> positive cross)
            cross = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
            if cross <= 0:
                continue
            # No other vertex inside this candidate ear.
            if any(p != i0 and p != i1 and p != i2 and
                   point_in_poly(poly[p], [a, b, c]) for p in idx):
                continue
            tris.append((i0, i1, i2))
            idx.pop(k)
            ear = True
            break
        if not ear:
            break
    if len(idx) == 3:
        tris.append((idx[0], idx[1], idx[2]))
    return tris


# ---------------------------------------------------------------------------
# Extrusion
# ---------------------------------------------------------------------------

def build_letter(letter):
    contours = glyph_contours(letter)
    # Normalise to unit height, centred at origin.
    allpts = [p for c in contours for p in c]
    xs = [p[0] for p in allpts]; ys = [p[1] for p in allpts]
    minx, maxx, miny, maxy = min(xs), max(xs), min(ys), max(ys)
    h = maxy - miny
    cx = 0.5 * (minx + maxx); cy = 0.5 * (miny + maxy)
    norm = [[((x - cx) / h, (y - cy) / h) for (x, y) in c] for c in contours]

    # Classify outer (largest |area|) vs holes; orient outer CCW, holes CW.
    areas = [abs(signed_area(c)) for c in norm]
    outer_i = max(range(len(norm)), key=lambda i: areas[i])
    outer = orient(norm[outer_i], ccw=True)
    holes = [orient(c, ccw=False) for i, c in enumerate(norm) if i != outer_i]

    simple = bridge_holes(outer, holes) if holes else outer
    simple = orient(simple, ccw=True)
    face_tris = ear_clip(simple)

    # Correctness gate: triangulated area == polygon area (outer - holes).
    poly_area = abs(signed_area(outer)) - sum(abs(signed_area(h)) for h in holes)
    tri_area = 0.0
    for (i, j, k) in face_tris:
        a, b, c = simple[i], simple[j], simple[k]
        tri_area += 0.5 * abs((b[0] - a[0]) * (c[1] - a[1])
                              - (b[1] - a[1]) * (c[0] - a[0]))
    assert abs(tri_area - poly_area) < 1e-3, \
        f'{letter}: triangulated area {tri_area:.5f} != polygon {poly_area:.5f}'

    verts = []; norms = []; tris = []
    zf = DEPTH / 2.0

    def add(v, n):
        verts.append(v); norms.append(n); return len(verts) - 1

    # Front (+z) and back (-z) faces.
    front_map = {}
    back_map = {}
    for vi, (x, y) in enumerate(simple):
        front_map[vi] = add((x, y, zf), (0.0, 0.0, 1.0))
    for vi, (x, y) in enumerate(simple):
        back_map[vi] = add((x, y, -zf), (0.0, 0.0, -1.0))
    for (i, j, k) in face_tris:
        tris.append((front_map[i], front_map[j], front_map[k]))
        tris.append((back_map[i], back_map[k], back_map[j]))   # reversed winding

    # Side walls from the ORIGINAL contours (outer + holes), each with its own
    # winding so the horizontal normal (dy,-dx) points out of the solid.
    for contour in [outer] + holes:
        m = len(contour)
        for e in range(m):
            x0, y0 = contour[e]
            x1, y1 = contour[(e + 1) % m]
            dx, dy = x1 - x0, y1 - y0
            ln = math.hypot(dx, dy)
            if ln < 1e-9:
                continue
            nx, ny = dy / ln, -dx / ln
            f0 = add((x0, y0, zf), (nx, ny, 0.0))
            f1 = add((x1, y1, zf), (nx, ny, 0.0))
            b0 = add((x0, y0, -zf), (nx, ny, 0.0))
            b1 = add((x1, y1, -zf), (nx, ny, 0.0))
            tris.append((f0, b0, f1))
            tris.append((f1, b0, b1))

    return (numpy.array(verts, dtype='float32'),
            numpy.array(norms, dtype='float32'),
            numpy.array(tris, dtype='int32'))


# ---------------------------------------------------------------------------
# Emit module
# ---------------------------------------------------------------------------

def fmt_array(name, arr):
    if arr.dtype == numpy.float32:
        rows = ['        [%s],' % ', '.join('%.4f' % v for v in row) for row in arr]
        dt = "'float32'"
    else:
        rows = ['        [%s],' % ', '.join('%d' % v for v in row) for row in arr]
        dt = "'int32'"
    return '%s = _np.array([\n%s\n    ], dtype=%s)\n' % (name, '\n'.join(rows), dt)


def main():
    out = []
    out.append('# @Author: Tristan Croll')
    out.append('# @License: Free for non-commercial use (see license.pdf)')
    out.append('# @Copyright: 2026 Tristan Croll')
    out.append('')
    out.append("'''")
    out.append('Precomputed extruded "R"/"S" letter meshes for ISOLDE chiral R/S markup.')
    out.append('')
    out.append('AUTO-GENERATED by ``tools/gen_rs_glyphs.py`` -- do not edit by hand.')
    out.append('Each letter is normalised to unit height, centred at the origin in the XY')
    out.append('plane, and extruded along z. ``letter_R_geometry()`` / ``letter_S_geometry()``')
    out.append('return ``(vertices, normals, triangles)`` matching ``Model.set_geometry``')
    out.append('(same contract as ``geometry.exclamation_mark``). No runtime font dependency.')
    out.append("'''")
    out.append('')
    out.append('import numpy as _np')
    out.append('')
    for letter, prefix in (('R', 'R'), ('S', 'S')):
        v, n, t = build_letter(letter)
        print(f'{letter}: {len(v)} verts, {len(t)} tris')
        out.append('# --- letter %s ---' % letter)
        out.append(fmt_array('_%s_VERTICES' % prefix, v))
        out.append(fmt_array('_%s_NORMALS' % prefix, n))
        out.append(fmt_array('_%s_TRIANGLES' % prefix, t))
    out.append('')
    out.append('def letter_R_geometry():')
    out.append("    '''(vertices, normals, triangles) for an extruded \"R\".'''")
    out.append('    return _R_VERTICES, _R_NORMALS, _R_TRIANGLES')
    out.append('')
    out.append('def letter_S_geometry():')
    out.append("    '''(vertices, normals, triangles) for an extruded \"S\".'''")
    out.append('    return _S_VERTICES, _S_NORMALS, _S_TRIANGLES')
    out.append('')
    with open(OUT_PATH, 'w') as f:
        f.write('\n'.join(out))
    print('wrote', OUT_PATH)


if __name__ == '__main__':
    main()
