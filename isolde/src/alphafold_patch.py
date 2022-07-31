# Temporary monkey-patch to ChimeraX 1.4 to handle v3 of the AlphaFold database
# REMOVE FOR 1.5

def alphafold_monkeypatch(session):
    from chimerax.alphafold.database import _alphafold_database_settings
    settings = _alphafold_database_settings(session)
    settings.database_version = 3

    def read_json_pae_matrix(path):
        '''Open AlphaFold database distance error PAE JSON file returning a numpy matrix.'''
        f = open(path, 'r')
        import json
        j = json.load(f)
        f.close()


        if isinstance(j, dict) and 'pae' in j:
            # ColabFold 1.3 produces a JSON file different from AlphaFold database.
            from numpy import array, float32
            pae = array(j['pae'], float32)
            return pae
        
        if not isinstance(j, list):
            from chimerax.core.errors import UserError
            raise UserError(f'JSON file "{path}" is not AlphaFold predicted aligned error data, expected a top level list')
        d = j[0]


        if not isinstance(d, dict):
            from chimerax.core.errors import UserError
            raise UserError(f'JSON file "{path}" is not AlphaFold predicted aligned error data, expected a top level list containing a dictionary')
            
        if 'residue1' in d and 'residue2' in d and 'distance' in d:
            # AlphaFold Database versions 1 and 2 use this format
            # Read PAE into numpy array
            from numpy import array, zeros, float32, int32
            r1 = array(d['residue1'], dtype=int32)
            r2 = array(d['residue2'], dtype=int32)
            ea = array(d['distance'], dtype=float32)
            # me = d['max_predicted_aligned_error']
            n = r1.max()
            pae = zeros((n,n), float32)
            pae[r1-1,r2-1] = ea
            return pae
            
        if 'predicted_aligned_error' in d:
            # AlphaFold Database version 3 uses this format.
            from numpy import array, float32
            pae = array(d['predicted_aligned_error'], dtype=float32)
            pae[pae==0]=0.2
            return pae
        
        keys = ', '.join(str(k) for k in d.keys())
        from chimerax.core.errors import UserError
        raise UserError(f'JSON file "{path}" is not AlphaFold predicted aligned error data, expected a dictionary with keys "predicted_aligned_error" or "residue1", "residue2" and "distance", got keys {keys}')

    from chimerax.alphafold import pae
    pae.read_json_pae_matrix = read_json_pae_matrix

    def alphafold_pae(session, structure = None, file = None, uniprot_id = None,
                    version = None, palette = None, range = None, plot = None,
                    color_domains = False, connect_max_pae = 5, cluster = 0.5, min_size = 10):
        '''Load AlphaFold predicted aligned error file and show plot or color domains.'''

        if uniprot_id:
            from chimerax.alphafold.database import alphafold_pae_url
            pae_url = alphafold_pae_url(session, uniprot_id, database_version=version)
            file_name = pae_url.split('/')[-1]
            from chimerax.core.fetch import fetch_file
            file = fetch_file(session, pae_url, 'AlphaFold PAE %s' % uniprot_id,
                            file_name, 'AlphaFold', error_status = False)
            
        if file:
            from chimerax.alphafold.pae import AlphaFoldPAE
            pae = AlphaFoldPAE(file, structure)
            if structure:
                if structure.num_residues != pae.matrix_size:
                    from chimerax.core.errors import UserError
                    raise UserError('Number of residues in structure "%s" is %d which does not match PAE matrix size %d.'
                                    % (str(structure), structure.num_residues, pae.matrix_size) +
                                    '\n\nThis can happen if the AlphaFold model has been trimmed to match an experimental structure, or if residues have been deleted.  The full-length AlphaFold model must be used to show predicted aligned error.')
                structure.alphafold_pae = pae
        elif structure is None:
            from chimerax.core.errors import UserError
            raise UserError('No structure or PAE file specified.')
        else:
            pae = getattr(structure, 'alphafold_pae', None)
            if pae is None:
                from chimerax.core.errors import UserError
                raise UserError('No predicted aligned error (PAE) data opened for structure #%s'
                                % structure.id_string)

        if plot is None:
            plot = not color_domains	# Plot by default if not coloring domains.
            
        if plot:
            from chimerax.core.colors import colormap_with_range
            from chimerax.alphafold.pae import AlphaFoldPAEPlot
            colormap = colormap_with_range(palette, range, default_colormap_name = 'pae',
                                        full_range = (0,30))
            p = getattr(structure, '_alphafold_pae_plot', None)
            if p is None or p.closed():
                p = AlphaFoldPAEPlot(session, 'AlphaFold Predicted Aligned Error', pae,
                                    colormap=colormap)
                if structure:
                    structure._alphafold_pae_plot = p
            else:
                p.display(True)
                if palette is not None or range is not None:
                    p.set_colormap(colormap)
                if divider_lines is not None:
                    p.show_chain_dividers(divider_lines)

        pae.set_default_domain_clustering(connect_max_pae, cluster)
        if color_domains:
            if structure is None:
                from chimerax.core.errors import UserError
                raise UserError('Must specify structure to color domains.')
            pae.color_domains(connect_max_pae, cluster, min_size)

    # -----------------------------------------------------------------------------
    #
    def register_alphafold_pae_command(logger):
        from chimerax.core.commands import CmdDesc, register, OpenFileNameArg, ColormapArg, ColormapRangeArg, BoolArg, FloatArg, IntArg
        from chimerax.atomic import AtomicStructureArg, UniProtIdArg
        desc = CmdDesc(
            optional = [('structure', AtomicStructureArg)],
            keyword = [('file', OpenFileNameArg),
                    ('uniprot_id', UniProtIdArg),
                    ('palette', ColormapArg),
                    ('range', ColormapRangeArg),
                    ('plot', BoolArg),
                    ('color_domains', BoolArg),
                    ('connect_max_pae', FloatArg),
                    ('cluster', FloatArg),
                    ('min_size', IntArg),
                    ('version', IntArg)],
            synopsis = 'Show AlphaFold predicted aligned error'
        )
        
        register('alphafold pae', desc, alphafold_pae, logger=logger)

    from chimerax.alphafold import pae
    pae.alphafold_pae = alphafold_pae
    pae.register_alphafold_pae_command = register_alphafold_pae_command

    def alphafold_fetch(session, uniprot_id, color_confidence=True,
                        align_to=None, trim=True, pae=False, ignore_cache=False,
                        add_to_session=True, version=None, in_file_history=True, **kw):

        uniprot_name = uniprot_id if '_' in uniprot_id else None
        from chimerax.alphafold.fetch import _parse_uniprot_id
        uniprot_id = _parse_uniprot_id(uniprot_id)
        from chimerax.alphafold import database
        url = database.alphafold_model_url(session, uniprot_id, version)
        file_name = url.split('/')[-1]
        
        from chimerax.core.fetch import fetch_file
        filename = fetch_file(session, url, 'AlphaFold %s' % uniprot_id, file_name, 'AlphaFold',
                            ignore_cache=ignore_cache, error_status = False)

        model_name = 'AlphaFold %s' % (uniprot_name or uniprot_id)
        models, status = session.open_command.open_data(filename, format = 'mmCIF',
                                                        name = model_name,
                                                        in_file_history = in_file_history,
                                                        **kw)
        from chimerax.alphafold.match import _set_alphafold_model_attributes
        _set_alphafold_model_attributes(models, uniprot_id, uniprot_name)

        if color_confidence:
            from chimerax.alphafold.fetch import _color_by_confidence
            for s in models:
                # Set initial style so confidence coloring is not replaced.
                s.apply_auto_styling()
                s._auto_style = False
                _color_by_confidence(s)

        if pae:
            trim = False	# Cannot associate PAE if structure is trimmed
            
        if align_to is not None:
            from chimerax.alphafold.fetch import _align_and_trim, _log_chain_info
            _align_and_trim(models, align_to, trim)
            _log_chain_info(models, align_to.name)
            
        if add_to_session:
            session.models.add(models)

        if pae:
            from chimerax.alphafold.pae import alphafold_pae
            alphafold_pae(session, structure = models[0], uniprot_id = uniprot_id)
            
        return models, status

    from chimerax.alphafold import fetch
    fetch.alphafold_fetch = alphafold_fetch

    def alphafold_sequence_search(sequences, min_length=20, local=False, log=None):
        '''
        Search all AlphaFold database sequences using blat.
        Return best match uniprot ids.
        '''
        from chimerax.alphafold.search import _search_sequences_local, _search_sequences_web, _plural
        useqs = list(set(seq for seq in sequences if len(seq) >= min_length))
        if len(useqs) == 0:
            return [None] * len(sequences)

        if log is not None:
            log.status('Searching AlphaFold database for %d sequence%s'
                    % (len(useqs), _plural(useqs)))

        if local:
            seq_uniprot_ids = _search_sequences_local(useqs)
        else:
            session.logger.warning('ChimeraX AlphaFold sequence search currently only searches v2 of the '
                'AlphaFold sequence database, covering the original release of 40 reference genomes. '
                'To search the complete database, go to https://www.ebi.ac.uk/Tools/sss/fasta/ in your web browser'
                'and tick (only) "AlphaFold DB" under "Structures" in the first window. Due to the '
                'size of the database, expect that to take a while.')
            seq_uniprot_ids = _search_sequences_web(useqs)
        seq_uids = [seq_uniprot_ids.get(seq) for seq in sequences]

        return seq_uids
    
    from chimerax.alphafold import search
    search.alphafold_sequence_search = alphafold_sequence_search
