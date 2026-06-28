# @Author: Tristan Croll
# @Date:   27-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 27-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Render endpoint for agent vision -- PARKED (UNEXPOSED) pending a Qt fix.

Goal: hand the agent (a) a render of what the human currently sees, plus (b) one
or more ORTHOGONAL views (different camera angles) that add spatial information
beyond the single on-screen view.

The implementation below is COMPLETE and correct; it is simply not registered on
the agent surface (see server.py ENABLE_RENDER / the commented MCP tool), because
of a **Qt 6.11 regression**: any offscreen image capture, after a tool window has
been floated, breaks the main window's resize-reshaping until a structural
relayout. This was bisected to the toolkit, not ISOLDE: it reproduces in *vanilla*
ChimeraX 1.13.dev (Qt 6.11) with `open <model>; save test.png width 300 height 200`
(after floating any panel), and does NOT occur in ChimeraX 1.12 (Qt 6.10). See
the write-up in scratch/chimerax_graphics_capture_report.md.

Re-enable (once Qt is fixed / a workaround lands): set ENABLE_RENDER = True in
server.py and uncomment the isolde_render tool in src/mcp/server.py, then verify
resize stays healthy after a capture.

Design notes for the live implementation:
- Capture runs INSIDE the redraw loop on the 'frame drawn' trigger (GL context
  current) and is bracketed with remember/restore_current_opengl_context.
- The user view is captured via session.main_view.image(); orthogonal views via a
  secondary View built on session.models.scene_root_model and initialised on the
  *same* OpenGLContext as the main view (so the scene's VAOs are valid -- a
  separate/offscreen context fails with GLError 1282, VAOs not being shared).
- `render` is registered raw (not thread_safe) so the HTTP thread waits on a Queue
  while the capture runs on the UI thread.
'''


def _agent_view(session):
    '''Cached secondary View for orthogonal-angle rendering, on the SAME GL
    context as the main view (shares the scene's GPU geometry incl. VAOs).'''
    main = session.main_view
    ctx = main.render.opengl_context
    av = getattr(session, '_isolde_agent_view', None)
    if av is not None and getattr(av, '_isolde_ctx', None) is ctx:
        return av
    from chimerax.graphics import View
    av = View(session.models.scene_root_model, window_size=(500, 400))
    av.initialize_rendering(ctx)
    av._isolde_ctx = ctx
    session._isolde_agent_view = av
    return av


def render_view(session, width=600, height=450, supersample=2,
                transparent_background=False, extra_angles=(90, 270), **_ignored):
    '''
    Render the human's current view plus optional orthogonal views, returned as
    base64 PNGs (MCP image content).

    Args:
        width, height: image size in pixels.
        supersample: anti-aliasing supersample factor (1 = off).
        transparent_background: render with an alpha channel.
        extra_angles: turntable angles (degrees, about the view-up axis through
            the centre of rotation) for additional orthogonal views; () for none.

    Returns {'images': [{label, image_base64, width, height, bytes}, ...]}.
    '''
    main = getattr(session, 'main_view', None)
    if main is None or getattr(main, 'render', None) is None:
        return {'error': 'no GL context (running headless --nogui?); use the GUI'}

    from queue import Queue, Empty
    from chimerax.core.triggerset import DEREGISTER
    q = Queue()

    def _encode(img, label):
        import io, base64
        buf = io.BytesIO()
        img.save(buf, format='PNG')
        data = buf.getvalue()
        return {'label': label, 'format': 'png', 'width': img.width,
                'height': img.height, 'bytes': len(data),
                'image_base64': base64.b64encode(data).decode('ascii')}

    def _capture(*_):
        from chimerax.graphics import (
            remember_current_opengl_context, restore_current_opengl_context)
        cc = remember_current_opengl_context()
        images = []
        try:
            # (a) the human's current view
            img = main.image(int(width), int(height), supersample=int(supersample),
                             transparent_background=bool(transparent_background))
            if img is not None:
                images.append(_encode(img, 'user'))
            # (b) orthogonal views via a secondary View on the SAME GL context
            if extra_angles:
                from chimerax.geometry import rotation
                av = _agent_view(session)
                try:
                    av.background_color = main.background_color
                    av.lighting = main.lighting
                    av.material = main.material
                except Exception:
                    pass
                cofr = main.center_of_rotation
                up = main.camera.position.axes()[1]   # view-up axis
                for ang in extra_angles:
                    av.camera.position = rotation(up, float(ang), cofr) * main.camera.position
                    im = av.image(int(width), int(height), supersample=int(supersample),
                                  transparent_background=bool(transparent_background))
                    if im is not None:
                        images.append(_encode(im, 'rot%d' % int(ang)))
            q.put({'images': images, 'count': len(images)})
        except Exception as e:
            import traceback
            q.put({'error': 'render failed: %s' % e, 'traceback': traceback.format_exc()})
        finally:
            restore_current_opengl_context(cc)
        return DEREGISTER

    def _arm():
        session.triggers.add_handler('frame drawn', _capture)
        main.redraw_needed = True

    session.ui.thread_safe(_arm)
    try:
        return q.get(timeout=15.0)
    except Empty:
        return {'error': 'render timed out waiting for a frame'}
