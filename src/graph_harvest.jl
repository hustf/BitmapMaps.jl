# Utilties concerning extracting data from the local and regional elevation graphs.
# - Instead of naive, recursive algorithms, we try to use DFSIterator.
# - Instead of repeated `inneighbor_labels` or `inneighbors`, reverse a copy of the graph
#   once and find descendants in it. 
# See `graph_utilties` for related functions.

"""
    harvest_summits_data_from_graph(gr)
    ---> Tuple{Vector{Float32}, Vector{Int64}, Vector{Tuple{Int64, Int64}}}

Output vectors: Elevation, summit prominence, utm position

gr is the regional graph gr. The local graphs will contain many 
leafs for which the prominence could not be determined, while
those have been pruned in the regional graph.
"""
function harvest_summits_data_from_graph(gr)
    vutm = leaf_utms(gr)
    vprom = prominence(gr, vutm)
    vz = map(utm -> round(Int64, gr[utm]), vutm)
    vz, vprom, vutm
end

leaf_utms(g) = [label_for(g, v) for v in vertices(g) if outdegree(g, v) == 0]

root_utms(g) = [label_for(g, v) for v in vertices(g) if indegree(g, v) == 0]

"""
    prominence(g, vutm)
    ---> Vector{Float32}

Calculate promincence for a vector of utm positions, which are labels in the graph g.
"""
function prominence(g, vutm)
    # Reversing the graph may speed up things if we can do it in bulk.
    # We would otherwise use 'inneigbors'. Now we can visit the same vertices,
    # but using 'outneighbors'.
    g_rev = reverse(g)
    f_down_iterator = let g_rev = g_rev
        (utm) -> DFSIterator(g_rev, code_for(g_rev, utm))
    end
    f_up_iterator = let g = g
        (utm) -> DFSIterator(g, code_for(g, utm))
    end
    map(vutm) do utm_top
        z_lim = g_rev[utm_top]
        iterator_nodes_below = f_down_iterator(utm_top)
        _prominence(g_rev, iterator_nodes_below, g, f_up_iterator, z_lim)
    end
end
function _prominence(g_rev, iterator_nodes_below, g, f_up_iterator, z_lim)
    for vertex_no in iterator_nodes_below
        utm = label_for(g_rev, vertex_no)
        z = g_rev[utm]
        candidate_prominence = z_lim - z
        candidate_prominence <= 0 && continue
        # Ok, is the valid candidate prominence also an actual prominence?
        if has_taller_descendant(g, f_up_iterator, z_lim, utm)
            return candidate_prominence
        end
    end
    # We found none taller than z_lim. There should be just one such summit per regional graph
    return z_lim
end


"""
    has_certain_smaller_prominence(g, min_prominence, vutm)
    ---> Vector{Bool}

Give no guarantee that the prominence can't be smaller than min_prominence, 
but a returned 'true' means that the prominence in now way can be larger.
"""
function has_certain_smaller_prominence(g, min_prominence, vutm)
    # Reversing the graph may speed up things if we can do it in bulk.
    # We would otherwise use 'inneigbors'. Now we can visit the same vertices,
    # but using 'outneighbors'.
    g_rev = reverse(g)
    f_down_iterator = let g_rev = g_rev
        (utm) -> DFSIterator(g_rev, code_for(g_rev, utm))
    end
    f_up_iterator = let g = g
        (utm) -> DFSIterator(g, code_for(g, utm))
    end
    map(vutm) do utm_top
        z_lim = g_rev[utm_top]
        iterator_nodes_below = f_down_iterator(utm_top)
        _has_certain_smaller_prominence(g_rev, iterator_nodes_below, g, f_up_iterator, z_lim, min_prominence)
    end
end
function _has_certain_smaller_prominence(g_rev, iterator_nodes_below, g, f_up_iterator, z_lim, min_prominence)
    for vertex_no in iterator_nodes_below
        utm = label_for(g_rev, vertex_no)
        z = g_rev[utm]
        candidate_prominence = z_lim - z
        # A lot of `utm_top` vertices will be linked to reference nodes at the same level, in which case
        # `candidate_prominence` will be zero when utm signifies the reference node.
        # If that reference node is connected to a taller point (taller than z_lim), it is meaningful to say 
        # that `utm_top` has a prominence of zero. 
        # There are two reasons for leaf nodes being added in the first place:
        # 1) This point `utm_top` was identified as having the right curvature and being a local maxium.
        #    We would not want to prune such a point. 
        # 2) `utm_top` was added because it might connect to another sheet.
        #    We want to prune such points if they failed to connect. However, those could
        #    be pruned from the `merge_into!` function without checking prominence.
        # Hence, we design this function for case 1)
        candidate_prominence < 1 && continue
        candidate_prominence >= min_prominence && continue
        # Ok, is the valid candidate prominence also an actual prominence?
        if has_taller_descendant(g, f_up_iterator, z_lim, utm)
            return true
        end
    end
    return false
end
function has_taller_descendant(g, f_up_iterator, z_lim, utm)
    up_iterator = f_up_iterator(utm)
    for vertex_up_no in up_iterator
        utm_up = label_for(g, vertex_up_no)
        z_up = g[utm_up]
        if z_up > z_lim
            # Yes, we found a taller position
            return true
        end
    end
    return false
end

function tallest_descendant(g, avoid_utm, utm)
    avoid = code_for(g, avoid_utm)
    start = code_for(g, utm)
    maximum(v -> v == avoid ? -9223372036854775808f0 : g[label_for(g, v)], DFSIterator(g, start))
end


function lowest_descendant(g_rev, utm)
    start = code_for(g_rev, utm)
    minimum(v -> g_rev[label_for(g_rev, v)], DFSIterator(g_rev, start))
end