# @Author: Tristan Croll <tic20>
# @Date:   11-Jun-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 20-Jul-2020
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
    ccd_name = template_map.get(template_name, None)
    if ccd_name is not None:
        return ccd_name
    if template_name.startswith('GLYCAM'):
        from .glycam import glycam_suffix_to_ccd_name
        ccd_name = glycam_suffix_to_ccd_name.get(template_name[-2:], None)
        if ccd_name is not None:
            return ccd_name
    if template_name.startswith('MC'):
        return template_name.split('_')[1]
    return None
