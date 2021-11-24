
def write_real_space_refine_defaults(session, model, volume, resolution, 
        param_file_name=None, model_file_name=None,
        restrain_coordination_sites=False,
        restrain_positions=False,
        include_hydrogens=False,
        ncs_constraints=False
        ):
    from chimerax.clipper.maps.map_handler_base import XmapHandlerBase
    if isinstance(volume, XmapHandlerBase):
        from chimerax.core.errors import UserError
        raise UserError('write_real_space_refine_defaults() should not be used for crystallographic data!')
    if param_file_name is None:
        param_file_name = 'real_space_refine.eff'
    if restrain_coordination_sites:
        raise NotImplementedError('Coordination site restraints not yet implemented')
    map_file = volume.path
    from chimerax.core.commands import run
    if model_file_name is None:
        model_file_name = '{}_for_phenix.cif'.format(model.name)
    if not include_hydrogens:
        sel_mask = model.atoms.selecteds
        model.atoms.selecteds=False
        model.atoms[model.atoms.element_names != 'H'].selecteds=True
        sel_text = 'sel true'
    else:
        sel_text = ''
    run(session, 'save {} #{} {}'.format(model_file_name, model.id_string, sel_text))
    if not include_hydrogens:
        model.atoms.selecteds=sel_mask
    out_str = f'''

data_manager {{
    real_map_files = "{map_file}"
    model {{
        file = "{model_file_name}"
        type = electron
    }}
}}

resolution = {resolution}
rotamers {{
    fit = all outliers_or_poormap *outliers_and_poormap outliers poormap
    restraints {{
        enabled = False
    }}
}}
write_initial_geo_file = False
write_final_geo_file = False
write_all_states = False
write_pkl_stats = False
model_format = pdb *mmcif
mask_and_he_map = False
resolution_factor = 0.25
ncs_constraints = {ncs_constraints}
refine_ncs_operators = Auto
variable_rama_potential = False
weight = None
real_space_target_exclude_selection = "element H or element D"
show_statistics = True
show_per_residue = True
dry_run = False
random_seed = 0
nproc = 12
gradients_method = fd linear quadratic *tricubic
scattering_table = n_gaussian wk1995 it1992 *electron neutron
ignore_symmetry_conflicts = False
wrapping = Auto
absolute_length_tolerance = 0.01
absolute_angle_tolerance = 0.01

refinement {{
    run = *minimization_global rigid_body local_grid_search morphing \
        simulated_annealing *adp
    nqh_flips = False
    max_iterations = 100
    macro_cycles = 5
    target_bonds_rmsd = 0.01
    target_angles_rmsd = 1.0
    backbone_sample = None
    use_adp_restraints = True
    atom_radius = 3
    do_ss_ideal = False
    do_ccd = False
    rigid_body {{
    group = None
    }}    
}}

reference_model {{
    enabled = True
    file = None
    use_starting_model_as_reference = True
    sigma = 1.0
    limit = 15.0
    hydrogens = False
    main_chain = True
    side_chain = True
    fix_outliers = False
    strict_rotamer_matching = True
}}

pdb_interpretation {{
    reference_coordinate_restraints {{
        enabled = {restrain_positions}
        exclude_outliers = False
        selection = "all"
        sigma = 0.2
        limit = 1.0
        top_out = True
    }}
    automatic_linking {{
        link_metals = True
        link_carbohydrates = True
        link_ligands = True
    }}
    peptide_link {{
        ramachandran_restraints = False
        restrain_rama_outliers = False
        restrain_rama_allowed = False
    }}
    ramachandran_plot_restraints {{
        enabled = False
    }}   
}}
    '''
    with open(param_file_name, 'wt') as outfile:
        outfile.write(out_str)
    import os
    session.logger.info(
        f'A phenix.real_space_refine input file {param_file_name} with settings recommended for ISOLDE models has been written to '
        f'{os.getcwd()} along with a current snapshot of your model ({model_file_name}). '
        f'You can start a refinement job by running the following command in the working directory: \n'
        f'phenix.real_space_refine {param_file_name}.'
    )



