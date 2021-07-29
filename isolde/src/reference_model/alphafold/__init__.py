
def fetch_alphafold_pae(session, uniprot_id, ignore_cache=False):
    '''
    Fetch the predicted aligned error matrix for an AlphaFold prediction from the 
    EMBL AlphaFold server.
    '''
    from chimerax.core.fetch import fetch_file
    file_name = f'AF-{uniprot_id}-F1-predicted_aligned_error_v1.json'
    url = f'https://alphafold.ebi.ac.uk/files/{file_name}'

    filename = fetch_file(session, url, f'Alphafold {uniprot_id} PAE', file_name, 'AlphaFold-PAE',
        ignore_cache=ignore_cache)
    
    import json, numpy

    with open(filename, 'rt') as f:
        data = json.load(f)[0]
    
    r1, r2, d = data['residue1'],data['residue2'],data['distance']

    size = max(r1)

    matrix = numpy.empty((size,size))

    matrix.ravel()[:] = d

    return matrix

