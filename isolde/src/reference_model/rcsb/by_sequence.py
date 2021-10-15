def find_matching_sequences(sequence, identity_cutoff=0.3, evalue_cutoff=1, dates_before="2021-03-01"):
    import requests
    from json import JSONDecodeError
    query = f'''{{
  "query": {{
    "type": "group",
    "logical_operator": "and",
    "nodes": [
      {{
        "type": "terminal",
        "service": "sequence",
        "parameters": {{
        "evalue_cutoff": {evalue_cutoff},
        "identity_cutoff": {identity_cutoff},
        "target": "pdb_protein_sequence",
        "value": "{sequence}"
        }}
      }},
      {{
          "type": "terminal",
          "service": "text",
          "parameters": {{
              "operator": "less",
              "value": "{dates_before}",
              "attribute": "rcsb_accession_info.initial_release_date" 
          }}
      }}
    ]
  }},
  "request_options": {{
    "scoring_strategy": "sequence"
  }},
  "return_type": "polymer_entity"
}}'''
    result = requests.get('https://search.rcsb.org/rcsbsearch/v1/query', params={'json': query})
    try:
        rdict = result.json()['result_set']
    except KeyError:
        print('No results found!')
        return None
    except JSONDecodeError:
        print('Decode error!')
        return None
    return rdict

def get_wwpdb_matches_for_all_chains(model, identity_cutoff=0.3, evalue_cutoff=1):
    seen = set()
    matches = {}
    for c in model.chains:
        if c.characters in seen:
            continue
        results = find_matching_sequences(c.characters, identity_cutoff, evalue_cutoff)
        if results is not None:
            matches[c.chain_id] = [r['identifier'] for r in results]
        else:
            matches[c.chain_id] = []
        seen.add(c.characters)
    return matches

def count_matches_for_pdb_list(session, pdb_ids, identity_cutoff=0.3, evalue_cutoff=1, log_file = 'pdb_matches.log'):
    from chimerax.open_command.cmd import provider_open
    with open(log_file, 'wt') as log:
        log.write('PDB ID,Chain ID,Match count,Matching ids\n')
        for pdb_id in pdb_ids:
            try:
                models = provider_open(session, [pdb_id])
                m = models[0]
            except:
                print(f'Failed on {pdb_id}!')
                continue
            matches = get_wwpdb_matches_for_all_chains(m, identity_cutoff, evalue_cutoff)
            for c, matching in matches.items():
                num_matches = len(matching)
                log.write(f'{pdb_id},{c},{num_matches},{" ".join(matching)}\n')
                log.flush()
            session.models.close(models)


