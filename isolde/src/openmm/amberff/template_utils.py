# @Author: Tristan Croll <tic20>
# @Date:   24-Jul-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 01-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

def _load_template_maps():
    import json, os
    data_dir = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(data_dir,'template_map_unique.json'), 'rt') as f:
        template_map = json.load(f)
    return template_map

template_map = _load_template_maps()

def template_name_to_ccd_name(template_name):
    ccd_name = None
    description = None
    ccd_name = template_map.get(template_name, None)
    if ccd_name is not None:
        return ccd_name, None
    if template_name.startswith('GLYCAM'):
        res_name = template_name.split('_')[1]
        from .glycam import (glycam_suffix_to_ccd_name,
            anchor_name_to_ccd_template,
            glycam_prefix_to_description
            )
        ccd_name = glycam_suffix_to_ccd_name.get(res_name[-2:], None)
        if ccd_name is None:
            ccd_name = anchor_name_to_ccd_template.get(res_name, None)
            description = "glycosylated"
        else:
            description = glycam_prefix_to_description.get(res_name[0], None)
        return ccd_name, description
    if template_name.startswith('MC'):
        return template_name.split('_')[1], None
    if template_name.startswith('USER'):
        return template_name.split('_')[1], None
    return None, None

_template_atom_to_ccd_atom_cache = {}

def match_template_atoms_to_ccd_atoms(session, template, ccd_name = None, timeout=5):
    if ccd_name is None:
        ccd_name, _ = template_name_to_ccd_name(template.name)
    if ccd_name is None:
        raise TypeError('MD template {} has no match in the CCD!'.format(template.name))
    matches = _template_atom_to_ccd_atom_cache.get((template.name, ccd_name))
    if matches is not None:
        return matches
    from chimerax.atomic import mmcif
    ccd_template = mmcif.find_template_residue(session, ccd_name)
    if len(template.atoms) > 2:
        # Try by name first (typically very fast)
        from chimerax.isolde.graph import make_graph_from_residue_template
        tgraph = template.graph
        ccd_graph = make_graph_from_residue_template(ccd_template)
        ti, ci, timed_out = tgraph.maximum_common_subgraph(ccd_graph, big_first=True, timeout=timeout)
        if timed_out:
            ti_2, ci_2, to_2 = tgraph.maximum_common_subgraph(ccd_graph, timeout=timeout)
            if len(ti_2) > len(ti):
                ti = ti_2
                ci = ci_2
                if to_2:
                    warn_str = ('Graph matching of MD template {} to CCD template {}'
                        ' timed out. Match may not be optimal').format(
                            template.name, ccd_name
                        )
                    session.logger.warning(warn_str)
        _template_atom_to_ccd_atom_cache[(template.name, ccd_name)] = (ti, ci)
        return ti, ci
    import numpy
    ti = []
    ci = []
    for i, a in enumerate(template.atoms):
        ti.append(i)
        candidates = []
        for j, ca in enumerate(ccd_template.atoms):
            if j not in ci:
                if ca.element.number == a.element.atomic_number:
                    candidates.append((j, ca))
        if not len(candidates):
            raise TypeError('MD Template {} does not match CCD template {}!'.format(template.name, ccd_template.name))
        elif len(candidates)==1:
            ci.append(candidates[0][0])
        else:
            # Check to see if we can distinguish by name
            name_match = False
            for j, ca in candidates:
                if ca.name == a.name:
                    name_match = True
                    candidates.append(j)
                    break
            if not name_match:
                # Just pick the first
                ci.append(candidates[0][0])
    ti = numpy.array(ti)
    ci = numpy.array(ci)
    _template_atom_to_ccd_atom_cache[(template.name, ccd_name)] = (ti, ci)
    return ti, ci
