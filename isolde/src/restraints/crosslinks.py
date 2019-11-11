def parse_crosslinks_file(session, filename):
    link_defs = []
    logger = session.logger
    with open(filename, 'rt') as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            if len(line.strip()):
                split_line = line.strip().split()
                if len (split_line) != 5:
                    logger.warning('Mis-formatted line: {}'.format(line))
                    continue
                try:
                    link_defs.append((*split_line[:2], *[float(f) for f in split_line[2:]]))
                except:
                    logger.warning('Failed to parse line: {}'.format(line))
    return link_defs

def apply_links(session, link_defs, model, strength, well_half_width, alpha,
    compressed_color=[255,0,128,255], ideal_color=[0,255,255,255],
    stretched_color=[255,0,0,255]):
    from chimerax.isolde import session_extensions as sx
    adrm = sx.get_adaptive_distance_restraint_mgr(model, 'Crosslinks')
    adrm.set_colormap(compressed_color, ideal_color, stretched_color)
    from chimerax.core.commands import AtomSpecArg
    parser = AtomSpecArg('crosslink defs')
    for link in link_defs:
        arg1, arg2, d, ldc, udc = link
        atoms = []
        err = False
        for arg in (arg1, arg2):
            aspec = parser.parse(arg, session)[0]
            atom = aspec.evaluate(session, models=[model]).atoms
            if len(atom) != 1:
                session.logger.warning('AtomSpec {} evaluates to more than one atom! Ignoring this link. '.format(arg))
                err = True
                break
            atoms.append(atom[0])
        if not err:
            lower_distance = d-ldc
            upper_distance = d-udc
            mean_distance = (lower_distance + upper_distance)/2
            tolerance = upper_distance-mean_distance
            c = well_half_width * 2 * tolerance

            adr = adrm.add_restraint(*atoms)
            adr.target=mean_distance
            adr.tolerance=tolerance
            adr.c = c
            adr.kappa = strength*(c/10)**2
            adr.alpha = alpha
            adr.enabled = True
