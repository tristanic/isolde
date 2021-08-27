
_pae_cache = dict()

def fetch_alphafold_pae(session, uniprot_id, ignore_matrix_cache=False, ignore_download_cache=False):
    '''
    Fetch the predicted aligned error matrix for an AlphaFold prediction from the 
    EMBL AlphaFold server.
    '''
    global _pae_cache
    if not ignore_matrix_cache:
        pae = _pae_cache.get(uniprot_id, None)
        if pae is not None:
            return pae
    from chimerax.core.fetch import fetch_file
    file_name = f'AF-{uniprot_id}-F1-predicted_aligned_error_v1.json'
    url = f'https://alphafold.ebi.ac.uk/files/{file_name}'

    filename = fetch_file(session, url, f'Alphafold {uniprot_id} PAE', file_name, 'AlphaFold-PAE',
        ignore_cache=ignore_download_cache)
    
    pae = _pae_cache[uniprot_id] = parse_pae_file(filename)
    return pae

def parse_pae_file(pae_json_file):
    import json, numpy

    with open(pae_json_file, 'rt') as f:
        data = json.load(f)[0]
    
    r1, r2, d = data['residue1'],data['residue2'],data['distance']

    size = max(r1)

    matrix = numpy.empty((size,size))

    matrix.ravel()[:] = d

    return matrix


def alphafold_id(model):
    metadata = model.metadata
    entry_keys = metadata.get('entry', None)
    if entry_keys is None:
        return None
    entry_data = metadata.get('entry data', None)
    if entry_data is None:
        return None
    try:
        id_index = entry_keys.index('id')
    except ValueError:
        return None
    import re
    result = re.match('AF\-[A-Z0-9]*\-F1', entry_data[id_index-1])
    if result is None:
        return None
    return result[0]