def write_phenix_refine_defaults(session, model, xmapset,
        file_name=None,
        restrain_coordination_sites=False,
        include_hydrogens=False,
        num_processors=1,
        num_macrocycles=6,
        nqh_flips=True,
        scattering_type='xray',
        ):
    if restrain_coordination_sites:
        raise NotImplementedError('Coordination site restraints not yet implemented')
    if file_name is None:
        file_name = 'refine.eff'
    model_file_name = '{}_for_phenix.cif'.format(model.name)
    reflections_file_name = '{}_for_phenix.mtz'.format(model.name)

    if not include_hydrogens:
        sel_mask = model.atoms.selecteds
        model.atoms.selecteds=False
        model.atoms[model.atoms.element_names != 'H'].selecteds=True
        sel_text = 'sel true'
    else:
        sel_text = ''
    from chimerax.core.commands import run
    run(session, 'save {} #{} {}'.format(model_file_name, model.id_string, sel_text))
    if not include_hydrogens:
        model.atoms.selecteds=sel_mask

    if scattering_type=='xray':
        scattering_type='n_gaussian'
    run(session, 'save {} #{}'.format(reflections_file_name, xmapset.id_string))
    out_str = f'''

refinement {{
    input {{
        pdb {{
            file_name = "{model_file_name}"
        }}
        xray_data {{
            file_name = "{reflections_file_name}"
            labels = "FOBS,SIGFOBS"
            high_resolution = None
            low_resolution = None
            r_free_flags {{
                file_name = "{reflections_file_name}"
                label = "Free_Flags"
                test_flag_value = {xmapset.live_xmap_mgr.free_flag}
            }}
        }}
    }}
    output {{
        write_eff_file = True
        write_geo_file = False
        write_final_geo_file = False
        write_def_file = False
        write_model_cif_file = True
        write_reflection_cif_file = False
        export_final_f_model = False
        write_maps = False
        write_map_coefficients = True
        write_map_coefficients_only = False
        pickle_fmodel = False
        n_resolution_bins = None
    }}
    electron_density_maps {{
        map_coefficients {{
            map_type = "2mFo-DFc"
            mtz_label_amplitudes = "2FOFCWT"
            mtz_label_phases = "PH2FOFCWT"
            fill_missing_f_obs = True
        }}
        map_coefficients {{
            map_type = "2mFo-DFc"
            mtz_label_amplitudes = "2FOFCWT_no_fill"
            mtz_label_phases = "PH2FOFCWT_no_fill"
        }}
        map_coefficients {{
            map_type = "mFo-DFc"
            mtz_label_amplitudes = "FOFCWT"
            mtz_label_phases = "PHFOFCWT"
        }}
    }}
    refine {{
        strategy = *individual_sites *individual_adp *occupancies
    }}
    main {{
        nqh_flips = {nqh_flips}
        number_of_macro_cycles = {num_macrocycles}
        nproc = {num_processors}
        scattering_table = {scattering_type}
    }}
    hydrogens {{
        refine = individual *riding Auto
    }}
    geometry_restraints.edits.dihedral.periodicity = None
    target_weights {{
        optimize_xyz_weight = True
        optimize_adp_weight = True
    }}
    reference_model {{
        enabled = True
        use_starting_model_as_reference = True
    }}
}}
    '''
    with open(file_name, 'wt') as outfile:
        outfile.write(out_str)
    import os
    session.logger.info(
        f'A phenix.refine input file {file_name} with settings recommended for ISOLDE models '
        f'has been written to {os.getcwd()} along with a current snapshot of your model ({model_file_name}) '
        f'and reflections ({reflections_file_name}). Note that the reflections file only contains '
        f'amplitudes - if you have a MTZ file with intensities you may wish to edit {file_name} '
        f'to substitute that. You can start a refinement job by running the following command in the working directory: \n'
        f'phenix.refine {file_name}.'
    )
