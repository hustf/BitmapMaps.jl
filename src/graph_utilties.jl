# Utilties concerning building and restructuring local and regional elevation graphs.
# Note that the MaxTree is a local elevation graph (and a DAG), but it is also a dense matrix
# which is too large to process for a whole region.
#
# The regional graph is indexed with utm indices instead of image indices.
#
# Maxtree-only methods like `leaf_indices`, `parent_of_leaf_indices` and
# `core_family_dictionary`, `distinct_summit_indices` are defined in `summits_on_sheet`.

"""
    build_and_save_sheet_graph(ffna_graph, maxtree, z, vI, f_I_to_utm, border_sides, min_prominence)
"""
function build_and_save_sheet_graph(ffna_graph, maxtree, z, vI, f_I_to_utm, border_sides, min_prominence)
    # We're building a meta graph, i.e. vertices encode geographical position and elevation as metadata and
    # meta-labels.
    # We build the graph from a pre-existing dense "matrix-graph", the maxtree, but we will
    # only include the vertices needed for finding topographical summit prominence in a region of
    # multiple sheets.
    # Borders between sheets are treated differently, since vertices on borders will connect
    # the different sheet graphs later.
    internal_bbox = bbox_internal(z, f_I_to_utm)
    # The first sheet in a region will have borders on the North and East sides only.
    # Expanding the internal_bbox on the South and West will reduce the vertex density on
    # those borders, because the matree has no vertices on the moved borders.
    limiting_bbox = grow_bbox_on_sides_excepting(internal_bbox, border_sides)
    # Graph to add into
    gl = get_new_unsaved_graph()
    # Add vertices, ancestors and descendants, including edges between
    n_vertices_added = add_indices_and_direct_lineage!(gl, maxtree, z, vI, f_I_to_utm)
    #
    # Now clean what we added in order to reduce size.
    #
    reduce_and_prune!(gl, limiting_bbox, min_prominence)
    #
    # Further cleanup must wait until merged with regional graph.
    # Save the local graph (in MGFormat by default)
    savegraph(ffna_graph, gl)
    gl
end

"""
    grow_bbox_on_sides_excepting(bb::@NamedTuple{min_x::Int64, min_y::Int64, max_x::Int64, max_y::Int64}, border_sides)
    ---> NamedTuple

This is used for specifying which sides of out MetaGraph we want to protect from pruning and reduction.
'bb' would describe the maximum extent of the graph (which is labelled by geographical position of each vertex).
'border_sides' specifies which sides we don't need to protect.
"""
function grow_bbox_on_sides_excepting(bb::@NamedTuple{min_x::Int64, min_y::Int64, max_x::Int64, max_y::Int64}, border_sides)
    min_x = bb.min_x - (:w ∉ border_sides)
    max_x = bb.max_x + (:e ∉ border_sides)
    min_y = bb.min_y - (:s ∉ border_sides)
    max_y = bb.max_y + (:n ∉ border_sides)
    (;min_x, min_y, max_x, max_y)
end


"""
    reduce_and_prune!(g, limiting_bbox, min_prominence)

Simplifies graph 'g' while preserving more info where vertex labels
border on 'limiting_bbox'.

We have two operations, 'reduce' (r) and 'prune' (p). The optimum order
of operations is something like '"rpprrpprppr", depending on the graph details.

The algoritm limits consecutive 'r' to one, and consecutive 'p' to two,
and exits when the last 'r' and 'p' had no effect.
"""
function reduce_and_prune!(g, limiting_bbox, min_prominence)
    n0 = nv(g)
    n3 = -1 # Scoping
    while true
        n1 = nv(g)
        chain_reduce!(g, limiting_bbox)
        nprunes = 0
        while true
            nprunes += 1
            n2 = nv(g)
            prune_low_prominence!(g, limiting_bbox, min_prominence)
            n3 = nv(g)
            n2 == n3 && break
            nprunes > 1 && break
        end
        n1 == n3 && break
    end
    @debug "            Vertex count after reduce and prune: $(round(n3 / n0, digits = 2)) of original"
    g
end



"""
    add_indices_and_direct_lineage!(g::MetaGraph, maxtree::MaxTree, z, vI::Vector{CartesianIndex{2}}, f_I_to_utm)

Modifies the graph g, and returns the total number of vertices added.

Starting at each index in vI, we add direct ascendants and direct descendants. The relationships are taken from maxtree,
and we modify the directed graph g.

# Definitions

Child A, of parent B:
   - identified by its index `iA` in the matrices `parent_indices` and 'z'.
   - `parent_indices[iA]` == iB, where iB ≠ iA.

Parent B of child A:
   - identified by its index `iB` in the matrices `parent_indices` and 'z'.
   - `parent_indices[iA]` == iB, where iB ≠ iA.

Ascendant C of A:
   - identified by its index `iC` in the matrices `parent_indices` and 'z'.
   - parent or recursive parent of A. A MaxTree has no loops, hence no need to check
     that iC ≠ iA.

Descendant D of B:
   - identified by its index `iD` in the matrices `parent_indices` and 'z'.
   - Child or recursive child of B. A MaxTree has no loops, hence no need to check
     that iD ≠ iB.

(A member of) direct lineage E of B:
   - identified by its index `iE` in the matrices `parent_indices` and 'z'.
   - ascendant or descendant of B. There is a path B <-> E, either through parent-> child
     or child -> parent, with no not reversal of direction.
"""
function add_indices_and_direct_lineage!(g::MetaGraph, maxtree::MaxTree, z, vI::Vector{CartesianIndex{2}}, f_I_to_utm)
    count = 0
    R = CartesianIndices(axes(z))
    Rlin = LinearIndices(axes(z))
    # Maxtree uses linear indices, so let's use it too.
    f_i_to_utm = i -> f_I_to_utm(R[i])
    # We already have the table for parents
    parent_indices = maxtree.parentindices
    # We store a lookup table for children as a dictionary.
    pc_dic = _parent_child_dictionary(parent_indices)
    for i in Rlin[vI]
        this_index_count = add_index_and_direct_ascendants!(g, parent_indices, z, f_i_to_utm, i)
        count += this_index_count
        that_index_count = add_index_and_direct_descendants!(g, pc_dic, z, f_i_to_utm, i)
        count += that_index_count
    end
    count
end

"""
    _parent_child_dictionary(parent_indices)
    ---> Dict{Int64, Set{Int64}}

The argument must be the parent_indices in a MaxTree.
"""
function _parent_child_dictionary(parent_indices)
    pc_dic = Dict{Int64, Set{Int64}}()
    for (ic, ip) in enumerate(parent_indices)
        if haskey(pc_dic, ip)
            push!(pc_dic[ip], ic)
        else
            push!(pc_dic, ip => Set{Int64}(ic))
        end
    end
    pc_dic
end

"""
add_index_and_direct_ascendants!(g, parent_indices, z, f_i_to_utm, i::Int)
    --> Int64

Note that this breaks the '!' convention, to also return the mutated argument.

Recursive function, modifies g and returns the number of vertices added.
"""
function add_index_and_direct_ascendants!(g, parent_indices, z, f_i_to_utm, i::Int)
    count = 0
    utm = f_i_to_utm(i)
    @assert utm isa Tuple{Int64, Int64}
    if ! haskey(g, utm)
        # Add utm as a loose vertex in the graph. The vertex is labelelled
        # with its position as a Tuple{Int, Int}, and vertex data is its elevation
        g[utm] = z[i]
        count += 1
    end
    # The maxtree parent of this vertex is most likely on a rounded same elevation and has index
    ip = parent_indices[i]
    # If no self-parenting
    if i !== ip
        utmp = f_i_to_utm(ip)
        # If the parent vertex is not already in the graph
        if ! haskey(g, utmp)
            # Recurse
            count += add_index_and_direct_ascendants!(g, parent_indices, z, f_i_to_utm, ip)
        end
        # If the edge between is not in the graph
        if ! haskey(g, utmp, utm)
            if ! haskey(g, utm, utmp)
                # Add the edge from parent (lower or roughly equal elevation) to child (this vertex)
                add_edge!(g, utmp, utm)
            else
                throw("this occurs for some reason!")
            end
        end
    end
    count
end



"""
add_index_and_direct_descendants!(g, parent_indices, z, f_i_to_utm, i::Int)
    --> Int64

Recursive function, modifies g and returns the number of vertices added.
"""
function add_index_and_direct_descendants!(g, pc_dic, z, f_i_to_utm, i::Int)
    if ! haskey(pc_dic, i)
        # Peaks do not have children. Borders often have.
        return 0
    end
    count = 0
    utm = f_i_to_utm(i)
    @assert utm isa Tuple{Int64, Int64}
    if ! haskey(g, utm)
        # Add utm as a loose vertex in the graph. The vertex is labelelled
        # with its position as a Tuple{Int, Int}, and vertex data is its elevation
        g[utm] = z[i]
        count += 1
    end
    # The maxtree children of this vertex is most likely on a higher elevation
    set_ic = pc_dic[i]
    for ic in set_ic
        if i !== ic
            utmc = f_i_to_utm(ic)
            # If the child vertex is not already in the graph
            if ! haskey(g, utmc)
                # Recurse
                count += add_index_and_direct_descendants!(g, pc_dic, z, f_i_to_utm, ic)
            end
            # If the edge between is not in the graph
            if ! haskey(g, utm, utmc)
                if ! haskey(g, utmc, utm)
                    # Add the edge from parent (lower or roughly equal elevation) to child
                    add_edge!(g, utm, utmc)
                else
                    throw("this occurs, for some reason!")
                end
            end
        end
    end
    count
end




"""
    get_graph(ffna_graph)

Load (or create, save and return ) the regional graph.

Usage of this graph type:

haskey(gr, ((56345, 6949028))) checks if a vertex is registered for this position.
gr[(56345, 6949028)] returns its value e.g. 43.2f0

"""
function get_graph(ffna_graph)::ElevationGraph
    # Look for it on disk.
    if isfile(ffna_graph)
        loadgraph(ffna_graph, "Elevation_graph", MGFormat())
    else
        g = get_new_unsaved_graph()
        # Save the new graph (in MGFormat by default)
        savegraph(ffna_graph, g)
        g
    end
end

"""
    get_new_unsaved_graph()
    ---> MetaGraph

Create the type of graph we need in this project in a type-stable manner.
"""
function get_new_unsaved_graph()
    # Note, we supply a named function, in order to avoid this:
    # ┌ Warning: Attempting to store MetaGraphsNext.var"#11#13".
    # │    JLD2 only stores functions by name.
    # │  This may not be useful for anonymous functions.
    MetaGraph(DiGraph(),                       # Directed graph
        Tuple{Int64, Int64},                   # Vertex label type
        Float32,                               # Vertex metadata type
        Nothing,                               # Edge metadata type
        Set{@NamedTuple{min_x::Int64, min_y::Int64, max_x::Int64, max_y::Int64}}(), # Graph metadata, mutable container.
        identity,  # The default is otherwise an anonymous function
        nothing    # Default edge value
        )
end

function _is_on_border(utm, borders)
    utm[1] == borders[:min_x] ||
    utm[1] == borders[:max_x] ||
    utm[2] == borders[:min_y] ||
    utm[2] == borders[:max_y]
end


"""
    prune_low_prominence!(g, min_prominence)
    prune_low_prominence!(g, limiting_bbox, min_prominence)
    prune_low_prominence!(g, limiting_bbox, min_prominence, vutm)
    ---> typeof(gr)

Removes the specified (or not) vertices and edges, starting at leafs and
pruning to first branch.
"""
function prune_low_prominence!(g, min_prominence)
    bb = bbox_internal(g)
    limiting_bbox = grow_bbox_on_sides_excepting(bb, Symbol[])
    prune_low_prominence!(g, limiting_bbox, min_prominence)
end
function prune_low_prominence!(g, limiting_bbox, min_prominence)
    # All leaves, except those on edges of limiting_bbox
    leaves_utm = filter(utm -> !_is_on_border(utm, limiting_bbox), leaf_utms(g))
    prune_low_prominence!(g, limiting_bbox, min_prominence, leaves_utm)
end
function prune_low_prominence!(g, limiting_bbox, min_prominence, vutm)
    has_small_proms = has_certain_smaller_prominence(g, min_prominence, vutm)
    vutm_with_smaller_prominence = [utm for (utm, has_small_prom) in zip(vutm, has_small_proms) if has_small_prom]
    #@debug "            Found $(length(vutm_with_smaller_prominence)) leaves with prominence certainly < $min_prominence out of $(length(vutm)) leaves"
    # Prune down to first branch or until we meet the limiting box.
    prune_to_first_branch!(g, limiting_bbox, vutm_with_smaller_prominence)
    # So far, we have deleted edges only. Now get rid of isolated vertices.
    remove_isolated_vertices!(g)
end

"""
    prune_to_first_branch!(gr, limiting_bbox, vutm)
    prune_to_first_branch!(gr, limiting_bbox, vertex_no::Int64)
    prune_to_first_branch!(gr, limiting_bbox, utm::Tuple{Int64, Int64})
    prune_to_first_branch!(gr, vertex_no::Int64)
    prune_to_first_branch!(gr, utm::Tuple{Int64, Int64})


Just snip the edges. Remove free vertices in bulk later.
"""
function prune_to_first_branch!(gr, limiting_bbox, vutm)
    foreach(utm -> prune_to_first_branch!(gr, limiting_bbox, utm), vutm)
    gr
end
prune_to_first_branch!(gr, limiting_bbox, utm::Tuple{Int64, Int64}) = prune_to_first_branch!(gr, limiting_bbox, code_for(gr, utm))
function prune_to_first_branch!(gr, limiting_bbox, vertex_no::Int64)
    # Caller knows there are no children / outneighbors.
    # Where sheets have been connected in a regional graph, some vertices will have multiple
    # parents / inneigbors. The graph is still acyclic.
    #if indegree(gr, vertex_no) !== 1
    #    throw(ArgumentError("Vertex no $(vertex_no) at utm $(label_for(gr, vertex_no)) has $(indegree(gr, vertex_no)) parents (inneigbors)"))
    #end
    parent_no = first(inneighbors(gr, vertex_no))
    # Snip this edge, leave the vertex unconnected. Vertex numbers remain unchanged.
    rem_edge!(gr, parent_no, vertex_no) || throw(ErrorException("Failed to remove edge (parent_no, vertex_no) = $((parent_no, vertex_no)) at utm $(label_for(gr, vertex_no))."))
    # Check the parent's edges.
    # We have removed our edge to parent already.
    if outdegree(gr, parent_no) > 0
        # Parent has other children. This, vertex_no, was the last pruned, the vertices are removed in `remove_isolated_vertices`.
        return 1
    elseif indegree(gr, parent_no) == 0
        # Parent is a root with other child(ren).  We have recursed unto the bottom. The snipped vertices are removed in `remove_isolated_vertices`.
        return 1
    else
        if _is_on_border(label_for(gr, parent_no), limiting_bbox)
            return 1
        else
            # Dive deeper recursively.
            return 1 + prune_to_first_branch!(gr, limiting_bbox, parent_no)
        end
    end
end
function prune_to_first_branch!(gr, vutm::Vector{Tuple{Int64, Int64}})
    foreach(utm -> prune_to_first_branch!(gr, utm),  vutm)
    gr
end
prune_to_first_branch!(gr, utm::Tuple{Int64, Int64}) = prune_to_first_branch!(gr, code_for(gr, utm))
function prune_to_first_branch!(gr, vertex_no::Int64)
    # Caller knows there are no children / outneighbors.
    parent_no = first(inneighbors(gr, vertex_no))
    # Snip this edge, leave the vertex unconnected. Vertex numbers remain unchanged.
    rem_edge!(gr, parent_no, vertex_no) || throw(ErrorException("Failed to remove edge (parent_no, vertex_no) = $((parent_no, vertex_no)) at utm $(label_for(gr, vertex_no))."))
    # Check the parent's edges.
    # We have removed our edge to parent already.
    if outdegree(gr, parent_no) > 0
        # Parent has other children. This, vertex_no, was the last pruned, the vertices are removed in `remove_isolated_vertices`.
        return 1
    elseif indegree(gr, parent_no) == 0
        # Parent is a root with other child(ren).  We have recursed unto the bottom. The snipped vertices are removed in `remove_isolated_vertices`.
        return 1
    else
        # Dive deeper recursively.
        return 1 + prune_to_first_branch!(gr, parent_no)
    end
end



"""
    remove_isolated_vertices!(gr)
    ---> typeof(gr)

Just so.
"""
function remove_isolated_vertices!(gr)
    # We use the Graphs.jl interface, not MetaGraphsNext
    # Collect isolated vertices into a vector
    isolated_vertices = [v for v in vertices(gr)
                         if indegree(gr, v) == 0 && outdegree(gr, v) == 0]
    # Remove isolated vertices in reverse order to avoid index issues
    for v in reverse(isolated_vertices)
        lab = label_for(gr, v)
        delete!(gr, lab)
    end
    gr
end


"""
    chain_reduce!(g, limiting_bbox)

Remove vertices with one in and one out connection,
provided that the vertex and its two neighbors are not on a border.
"""
function chain_reduce!(g, limiting_bbox)
    # We'll remove superfluous vertices,
    # i.e. those with one neighbor in each direction.
    remove_links = Vector{Tuple{Int64, Int64}}()
    for v in reverse(vertices(g))
        if indegree(g, v) == 1 && outdegree(g, v) == 1
            utm = label_for(g, v)
            if ! _is_on_border(utm, limiting_bbox)
                utm1 = first(inneighbor_labels(g, utm))
                if ! _is_on_border(utm1, limiting_bbox)
                    utm2 = first(outneighbor_labels(g, utm))
                    if ! _is_on_border(utm2, limiting_bbox)
                        push!(remove_links, utm)
                    end
                end
            end
        end
    end
    # We're deleting the centre link from each triplet. The order of triplets is already optimum.
    for utm in remove_links
        # We cant rely on stored neighbors, since they may have been removed already.
        utmin = first(inneighbor_labels(g, utm))
        utmout = first(outneighbor_labels(g, utm))
        # Add shortcut over the cut vertex utm
        add_edge!(g, utmin, utmout)
        # Remove edges
        delete!(g, utmin, utm)
        delete!(g, utm, utmout)
        # Remove vertex
        delete!(g, utm)
    end
end


"""
    merge_into!(gr, gl, cell_iter, f_I_to_utm, border_sides)

gr is the regional (larger) elevation graph
gl is the local graph

border_sides indicates which sides of the local graph 'gl' have borders with
any sheet.

, gl. E.g. pairs_of_border_neighbors(gr, cell_iter, f_I_to_utm, border_sides)

Connects edges where relevant, enabling longer downhill and uphill paths.

Of course, labels from different sheets are necessarily unique,
describing geographical position. For the same reason, edges are unique.
"""
function merge_into!(gr, gl, cell_iter, f_I_to_utm, border_sides)
    _merge_into!(gr, gl)
    # Identify (likely still unconnected) pairs of neighbor coordinates.
    # The first coordinate in each pair is from the local graph, the second
    # existed previously
    @debug "            Find connecting pairs"
    bp = pairs_of_border_neighbors(gr, cell_iter, f_I_to_utm, border_sides)
    if isempty(bp)
        @debug "            No connecting pairs, no connections made while merging"
    end
    # Pre-determine the new edges directions,  based on the still un-connected graphs.
    # We bypass destination vertices that are leaf nodes, and instead
    # connect to it's reference node. That way, we can follow paths uphill across the
    # border while iteratively going source -> destination.
    #
    @debug "            Find source => destination for $(length(bp)) pairs"
    bp_directed = directed_pairs_of_border_neighbors(gl, gr, bp)
    # Now change the regional graph gr by adding edges across the border(sides).
    @debug "            Making $(length(bp_directed)) edges across border. $(length(bp) - length(bp_directed)) connections discarded"
    for (s, d) in bp_directed
        make_connection!(gr, s, d)
    end
    # Clean up: Remove leaf nodes still remaining on the border. Most of these occur because a reference node was targeted
    # instead of the vertex on the border.
    # It would be difficult to determine the prominence of leaves on the border at this point.
    # However, the likelihood of actual summits with prominence > min_prominence is miniscule.
    # We simply prune all branches with leaves on the border.
    while true
        discard_leaves = [utm for pair in bp for utm in pair if length(outneighbor_labels(gr, utm)) == 0 && length(inneighbor_labels(gr, utm)) == 1]
        isempty(discard_leaves) && break
        @debug "            Pruning $(length(discard_leaves)) branches ending on the border"
        prune_to_first_branch!(gr, discard_leaves)
    end
    # So far, we have deleted edges only. Now get rid of isolated vertices.
    remove_isolated_vertices!(gr)
    # Check that this is still a DAG.
    @debug "            Checking for loops"
    if is_cyclic(gr)
        throw(ErrorException("Inadvertently created a loop! "))
    end
    # Add this sheet to the regional graph's metadata
    internal_bbox = bbox_internal(cell_iter, f_I_to_utm)
    push!(gr[], internal_bbox)
    #
    gr
end

function make_connection!(gr, s, d)
    if ! add_edge!(gr, s, d)
        @show haskey(gr, s, d)
        throw(ErrorException("Failed to add edge s =$s , d = $d. Consider deleting the regional graph *.z file!"))
    end
    gr
end

function destination_proxy(gr, d::Tuple{Int64, Int64})
    if outdegree(gr, code_for(gr, d)) == 0
        if indegree(gr, code_for(gr, d)) == 1
            refnode_cand = first(inneighbor_labels(gr, d))
            if outdegree(gr, code_for(gr, refnode_cand)) > 0
                # This is the one we ought to connect to instead
                if round(Int64, gr[d]) == round(Int64, gr[refnode_cand])
                    return refnode_cand
                end
            else
                throw(ErrorException("Unexpected. Consider d $d."))
            end
        else
            throw(ErrorException("Unexpected, consider d = $d."))
        end
    end
    return d
end


"""
    _merge_into!(gr, gl)

Merges vertices, edges and graph metadata from gl into gr.
"""
function _merge_into!(gr, gl)
    # Forget about vertex numbers (codes) in gl, collect data and labels
    dic = Dict{Tuple{Int64, Int64}, Float32}()
    for lab in labels(gl)
        push!(dic, lab => gl[lab])
    end
    # Collect edges in terms of labels
    set = Set{Tuple{Tuple{Int64, Int64}, Tuple{Int64, Int64}}}()
    for e in edge_labels(gl)
        push!(set, e)
    end
    # Incrementally add vertices
    for (lab, val) in dic
        gr[lab] = val
    end
    # Incrementally add edge
    for (s, d) in set
        add_edge!(gr, s, d)
    end
    # Graph metadata
    for data in gl[]
        push!(gr[], data)
    end
    gr
end





"""
    pairs_of_border_neighbors(gr, cell_iter, f_I_to_utm, border_sides)
    pairs_of_border_neighbors(gr, cell_iter, f_I_to_utm, border_side::Symbol)
    --> Vector{Pair{Tuple{Int64, Int64}, Tuple{Int64, Int64}}}

Return pairs of coordinates which has vertices in graph gr.

The first coordinate in a pair is internal to R_utm, the second is external to it on border_side.

Border_sides are directions :n, :w, etc.
"""
function pairs_of_border_neighbors(gr, cell_iter, f_I_to_utm, border_sides)
    # Output container
    pairs = Vector{Pair{Tuple{Int64, Int64}, Tuple{Int64, Int64}}}()
    # Now connect pairs where possible
    foreach(border_sides) do border_side
        # Find the neighbouring, still unconnected vertices at this border side
        append!(pairs, pairs_of_border_neighbors(gr, cell_iter, f_I_to_utm, border_side))
    end
    pairs
end
function pairs_of_border_neighbors(gr, cell_iter, f_I_to_utm, border_side::Symbol)
    # Prepare
    cell2utm = cell_to_utm_factor(f_I_to_utm)
    Δ = neighbor_delta_utm(cell2utm, border_side)
    # Border side refers to the geographical area defined by f_I_to_utm.(cell_iter).
    fullside = indices_of_border(cell_iter, border_side)
    internalside = fullside[2:end - 1]
    border_utms = f_I_to_utm.(internalside)
    # Output container
    pairs = Vector{Pair{Tuple{Int64, Int64}, Tuple{Int64, Int64}}}()
    # Append pairs of utm coordinates, provided both exist in gr
    for internal_neighbor in border_utms
        if ! haskey(gr, internal_neighbor)
            throw(ArgumentError("Internal border on :$(border_side) is missing, at utm $(internal_neighbor)"))
        end
        external_neighbor = internal_neighbor .+  Δ
        # The neigbouring coordinate may or may not have been added to the graph already.
        if haskey(gr, external_neighbor)
            push!(pairs, internal_neighbor  => external_neighbor)
        end
    end
    pairs
end

"""
    neighbor_delta_utm(cell2utm, border_side)
    ---> Tuple{Int64, Int64}
"""
function neighbor_delta_utm(cell2utm, border_side)
    @assert cell2utm > 0
    if :n == border_side
        (0, cell2utm)
    elseif :s == border_side
        (0, -cell2utm)
    elseif :w == border_side
        (-cell2utm, 0)
    elseif :e == border_side
        (cell2utm, 0)
    else
        throw(ArgumentError("border_side is $border_side, not :n, :s, :e:, :w"))
    end
end


tallest_elevation_above(g, utm) = tallest_descendant(g, utm, destination_proxy(g, utm))



"""
    directed_pairs_of_border_neighbors(g1, g2, bp)
    ---> Vector{Pair{Tuple{Int64, Int64}, Tuple{Int64, Int64}}}

g1 and g2 are neighbouring elevation graphs.
bp is a vector with (utm1 => utm2)
    where utm1 is a border position in g1
          utm2 is the opposite border position in g2


Order the connections' endpoints as source => destination, depending on
which position is at a higher elevation, and other properites of the graphs.

If the destination is a leaf in a connected component, it is
replaced with its reference node. See maxtree definitions.

If a connection is symmetric, loosely speaking, this is filtered out.
Otherwise, the result and bp will be of the same length.
"""
function directed_pairs_of_border_neighbors(g1, g2, bp)
    # Prepare
    #
    # We need to move towards ascdendants many times,
    # which is slower than towards descendants.
    # An equivalent operation is to reverse the graph, then
    # walk towards ascendants. We reverse the graph once
    # and capture the reversed graphs in these closures:
    lowest_elevation_below_1 = let g1r = reverse(g1)
        (utm) -> lowest_descendant(g1r, utm)
    end
    lowest_elevation_below_2 = let g2r = reverse(g2)
        (utm) -> lowest_descendant(g2r, utm)
    end
    # CONSIDER
    # This is terribly slow at times. It might be easier to build
    # a 'cumulative' graph, where each vertex data contained it's own
    # lowest descendant value, and then use that for lookup.
    _directed_pairs_of_border_neighbors = let g1 = g1, g2 = g2
        (utm_pair) -> let
            utm1, utm2 = utm_pair
            z_pair = round(Int64, g1[utm1]) => round(Int64, g2[utm2])
            if z_pair == (0 => 0)
                return (nothing, nothing) => (nothing, nothing)
            else
                tallest_above_pair = round(Int64, tallest_elevation_above(g1, utm1)) => round(Int64, tallest_elevation_above(g2, utm2))
                lowest_below_pair = round(Int64, lowest_elevation_below_1(utm1)) => round(Int64, lowest_elevation_below_2(utm2))
                s, d = source_destination(utm_pair, z_pair, tallest_above_pair, lowest_below_pair)
                if d == utm1
                    dp = destination_proxy(g1, utm1)
                elseif d == utm2
                    dp = destination_proxy(g2, utm2)
                else
                    @assert s == (nothing, nothing) && d == (nothing, nothing)
                    dp = d
                end
                return s => dp
            end
        end
    end
    #
    dbp = map(_directed_pairs_of_border_neighbors, bp)
    filter(p -> p !== ((nothing, nothing) => (nothing, nothing)), dbp)
end


"""
    source_destination(utm_pair, z_pair, tallest_above_pair, lowest_below_pair)
    -> Pair{Tuple{Int64, Int64}, Tuple{Int64, Int64}}

Returns `(source_utm => destination_utm)`, an edge which will connect
the graphs without creating loops. Don't use it directly as an edge! Replace
`destination_utm` with its reference node first. This way, you can walk downhil through
ascendants only.

`source_utm` will be utm1 or utm2
`destination_utm` will be the remaining utm1 or utm2.
"""
function source_destination(utm_pair, z_pair, tallest_above_pair, lowest_below_pair)
    utm1, utm2 = utm_pair
    # Which coordinate is on a higher elevation?
    # The maxtree is made from rounded elevations, so we use that here too.
    if >(z_pair...)
        # The higher, one, will be the child, i.e. the destination vertex
        d = utm1
        s = utm2
    elseif <(z_pair...)
        s = utm1
        d = utm2
    else
        # s and d are on the same elevation.
        # Cases considering highest elevation above and lowest elevation below:
        # 1) utm1 or its reference node has tallest elevation above
        #    => s = utm2, d = utm1
        # 2) utm2 or its reference node has tallest elevation above
        #    => s = utm1, d = utm2
        # 3) both (reference) nodes have the same tallest elevation above
        #   3a) utm1 has lowest elevation below
        #    => s = utm1, d = utm2
        #   3b) utm2 has lowest elevation below
        #    => s = utm2, d = utm1
        #   3c) both have the same lowest elevation below
        #    => Don't make the connection.
        #
        # Note alternative approach: Link directly to the tallest descendant or lowest ascendant, not
        # allowing border crossing. This should work well if we only had two regions, but we're stitching
        # others later.
        if >(tallest_above_pair...) # 1)
            s = utm2
            d = utm1
        elseif <(tallest_above_pair...) # 2)
            s = utm1
            d = utm2
        else # 3)
            if <(lowest_below_pair...) # 3a)
                s = utm1
                d = utm2
            elseif >(lowest_below_pair...) # 3b)
                s = utm2
                d = utm1
            else # 3c)
                # Don't connect
                s = (nothing, nothing)
                d = (nothing, nothing)
            end
        end
    end
    return s => d
end
