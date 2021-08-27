
def cluster_into_domains(session, model, pae_matrix, pae_power=1, pae_cutoff=5, adjust_weights_for_distance=True, distance_power=1, graph_resolution=1, color_by_cluster=True):
    import networkx as nx
    import numpy
    weights = 1/pae_matrix**pae_power
    if adjust_weights_for_distance:
        distances = distance_matrix(model, num_residues = pae_matrix.shape[0])
        weights = weights/distances**distance_power
    # Within the NetworkX greedy_modularity_communities() method the weight is used to define a term
    # that becomes one element of a dict key. Therefore the use of floating point is inadvisable (can 
    # lead to random `KeyError`s). So, we convert the weights to integers, first multiplying by a 
    # sufficiently large number to make sure everything is left of the decimal place.
    # Note that as of 19 August 2021 `KeyError`s can still happen - this has been reported to the
    # NetworkX developers (https://github.com/networkx/networkx/issues/4992) and possible fixes are being 
    # explored.
    #weights = (weights * 1e6).astype(numpy.int)

    g = nx.Graph()
    size = weights.shape[0]
    g.add_nodes_from(range(size))
    for i in range(size):
        for j in range(i+1, size):
            pae = pae_matrix[i,j]
            if pae < pae_cutoff:
                weight = weights[i,j]
                if weight > 0:
                    g.add_edge(i, j, weight=weight)

    from networkx.algorithms import community

    clusters = community.greedy_modularity_communities(g, weight='weight', resolution=graph_resolution)
    residue_clusters = []
    for i, c in enumerate(clusters):
        residues = model.residues[numpy.in1d(model.residues.numbers,[r+1 for r in c])]
        residue_clusters.append(residues)
        for r in residues:
            r.isolde_domain = i
    
    if color_by_cluster:
        from chimerax.core.commands import run
        run(session, f'color byattribute r:isolde_domain #{model.id_string} target ra palette Paired-12')

    return residue_clusters

def distance_matrix(m, num_residues = None):
    if num_residues is None:
        num_residues = max(m.residues.numbers)
    from chimerax.geometry import distance
    import numpy
    size = num_residues
    distances = numpy.ones([size,size])
    cas = m.atoms[m.atoms.names=='CA']
    for i, ca1 in enumerate(cas):
        for ca2 in cas[i+1:]:
            dist = distance(ca1.coord, ca2.coord)
            r1num = ca1.residue.number
            r2num = ca2.residue.number
            distances[r1num-1,r2num-1] = dist
            distances[r2num-1,r1num-1] = dist
    return distances


def domains_from_pae_matrix(pae_matrix, pae_power=1, pae_cutoff=5, graph_resolution=1):
    '''
    Takes a predicted aligned error (PAE) matrix representing the predicted error in distances between each 
    pair of residues in a model, and uses a graph-based community clustering algorithm to partition the model
    into approximately rigid groups.

    Arguments:

        * pae_matrix: a (n_residues x n_residues) numpy array. Diagonal elements should be set to some non-zero
          value to avoid divide-by-zero warnings
        * pae_power (optional, default=1): each edge in the graph will be weighted proportional to (1/pae**pae_power)
        * pae_cutoff (optional, default=5): graph edges will only be created for residue pairs with pae<pae_cutoff
        * graph_resolution (optional, default=1): regulates how aggressively the clustering algorithm is. Smaller values
          lead to larger clusters. Value should be larger than zero, and values larger than 5 are unlikely to be useful.

    Returns: a series of lists, where each list contains the indices of residues belonging to one cluster.
    '''
    import networkx as nx
    import numpy
    weights = 1/pae_matrix**pae_power
    # Within the NetworkX greedy_modularity_communities() method the weight is used to define a term
    # that becomes one element of a dict key. Therefore the use of floating point is inadvisable (can 
    # lead to random `KeyError`s). So, we convert the weights to integers, first multiplying by a 
    # sufficiently large number to make sure everything is left of the decimal place.
    # Note that as of 19 August 2021 `KeyError`s can still happen - this has been reported to the
    # NetworkX developers (https://github.com/networkx/networkx/issues/4992) and possible fixes are being 
    # explored.
    weights = (weights * 1e6).astype(numpy.int)

    g = nx.Graph()
    size = weights.shape[0]
    g.add_nodes_from(range(size))
    for i in range(size):
        for j in range(i+1, size):
            pae = pae_matrix[i,j]
            if pae < pae_cutoff:
                weight = weights[i,j]
                if weight > 0:
                    g.add_edge(i, j, weight=weight)

    from networkx.algorithms import community

    clusters = community.greedy_modularity_communities(g, weight='weight', resolution=graph_resolution)
    return clusters
