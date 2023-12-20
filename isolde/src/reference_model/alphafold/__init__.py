
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

def b_factors_look_like_plddt(m):
    import numpy
    if m.atoms.bfactors.max() > 100:
        return (False, 'out of range')
    from chimerax.atomic import Residue, Atoms
    coil = m.residues[numpy.logical_and(m.residues.polymer_types==Residue.PT_AMINO, m.residues.ss_types==Residue.SS_COIL)]
    coil_median = numpy.median(Atoms(coil.principal_atoms).bfactors)
    structured = m.residues[numpy.logical_and(m.residues.polymer_types==Residue.PT_AMINO, m.residues.ss_types!=Residue.SS_COIL)]
    structured_median = numpy.median(Atoms(structured.principal_atoms).bfactors)
    if structured_median < coil_median:
        return (False, 'inverted')
    for r in m.residues:
        if not numpy.allclose(r.atoms.bfactors, r.atoms.bfactors.mean()):
            return (False, 'different within residues')
    return (True, None)

plddt_failure_reasons = {
    'out of range': 'Some b-factors are over 100, but pLDDT values should be provided on a 0-100 or 0-1 scale.',
    'inverted': 'High pLDDT values correspond to high confidence, but the median value for unstructured residues is higher than the median value for structured ones.',
    'different within residues': 'For models with pLDDT scores in the B-factor column, all atoms in a residue typically have the same score assigned. That is not true for this model.'
}

def convert_plddt_to_b_factors(m):
    from chimerax.core.errors import UserError
    looks_like, reason = b_factors_look_like_plddt(m)
    if not looks_like:
        raise UserError('The values in the B-factor column of this model do not look like pLDDT values. '
            f'{plddt_failure_reasons[reason]}')
    bvals = m.atoms.bfactors
    # Normalise to (0..1)
    if max(bvals) > 1:
        bvals /= 100
    import numpy
    from math import pi
    # From the supplementary info of https://doi.org/10.1126/science.abj8754
    err = 1.5*numpy.exp(4*(0.7-bvals))
    m.atoms.bfactors = 8*pi**2*err*2/3
    
def pae_from_pkl_to_json(pae_file, precision=3):
    '''
    The standalone version of AlphaFold only writes a pickle file with results, rather than the .json file provided by the AlphaFold DB. I am reluctant
    to support opening pickle files from the GUI since malicious pickle files can run arbitrary code on opening - but this method can be used to extract 
    the PAE to .json where needed.
    '''
    import pickle, json, os
    with open(pae_file, 'rb') as infile:
        pkl_dict = pickle.load(infile)
    out_name = os.path.splitext(pae_file)[0]+"_pae.json"
    pae = pkl_dict['predicted_aligned_error']
    rounded_pae = [[round(r, 3) for r in row] for row in pae]
    with open(out_name, 'wt') as outfile:
        json.dump([{'predicted_aligned_error': rounded_pae}], outfile)
