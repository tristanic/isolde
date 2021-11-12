def dbscan(points, epsilon=3, min_points=6):
    '''
    Simple implementation of the DBSCAN algorithm for clustering points in 3D
    space. 
    '''
    import numpy
    from chimerax.geometry import find_close_points
    core = set()
    neighbors = []
    for i, p in enumerate(points):
        _, close_i = find_close_points([points[i]], points, max_distance=epsilon)
        neighbors.append(set(close_i))
        if len(close_i) > min_points:
            core.add(i)
    
    remainder = set(range(len(points)))
    clusters = [] # final clusters
    # Assign core groups (those with at least min_points points) to clusters
    while len(core) > 0:
        remainder_old = remainder.copy()
        idx = core.pop()
        qp = [idx]
        remainder.remove(idx)
        while len(qp) > 0:
            q = qp.pop(0)
            nq = neighbors[q]
            if len(nq) >= min_points:
                delta = nq.intersection(remainder)
                qp.extend(delta)
                remainder.difference_update(delta)
        new_cluster = remainder_old.difference(remainder)
        clusters.append(new_cluster)
        core.difference_update(new_cluster)
    
    
    clusters = list(sorted(clusters, key=lambda c: len(c), reverse=True))

    # Add remaining points to clusters
    for cluster in clusters:
        close = [1]
        while len(close):
            remaining = numpy.array(list(remainder))
            rpoints = points[remaining]
            idx_map = {i: idx for i, idx in enumerate(remaining)}
            cpoints = points[numpy.array(list(cluster))]
            _, close = find_close_points(cpoints, rpoints, epsilon)
            if len(close):
                ci = set([idx_map[i] for i in close])
                cluster.update(ci)
                remainder.difference_update(ci)

    clusters = list(sorted(clusters, key=lambda c: len(c), reverse=True))
    clusters = [numpy.array(list(c)) for c in clusters]
    noise = numpy.array(list(remainder))
    return clusters, noise
