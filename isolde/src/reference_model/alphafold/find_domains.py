
def cluster_into_domains(model, pae_matrix, distance_cutoff=5, graph_resolution=1, color_by_cluster=True):
    import networkx as nx
    import numpy
    weights = (1/pae_matrix*10000).astype(numpy.int)
    g = nx.Graph()
    size = weights.shape[0]
    g.add_nodes_from(range(size))
    for i in range(size):
        for j in range(i, size):
            pae = pae_matrix[i,j]
            if pae < distance_cutoff:
                g.add_edge(i, j, weight=weights[i,j])

    from networkx.algorithms import community

    clusters = community.greedy_modularity_communities(g, weight='weight', resolution=graph_resolution)
    residue_clusters = []
    for i, c in enumerate(clusters):
        residues = model.residues[numpy.array(list(c))]
        residue_clusters.append(residues)
        for r in residues:
            r.isolde_domain = i
    
    if color_by_cluster:
        from chimerax.core.commands import run
        run(session, f'color byattribute r:isolde_domain #{model.id_string} target ra palette rainbow')

    return residue_clusters

