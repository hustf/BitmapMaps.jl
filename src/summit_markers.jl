# This file concerns finding summit prominence.
# The prominence of some summits depends
# on geography which lay outside the limits of one sheet.
# We can't fit every sheet in memory at once, so 
# a boundary interaction solution should be found after al the initial single sheet
# calculations are finished. 
# This file does not concern that interaction,
# other than it prepares and saves an intermediate step, namely
# the 'maximum elevation above' (mea) table.

"""
    summit_markers(sb::SheetBuilder)
    summit_markers(fofo, cell_iter, cell2utm, f_I_to_utm, prom_levels, symbols, symbol_sizes)
    --> Bool

- Identify 'prominent' and 'obscure' summits.
- Mark symbols in image MARKERS_FNAM. 
- Summarize in .csv file for potential name tagging: SUMMITS_FNAM.
- Also writes the intermediate matrix 'maximum elevation above' to MAX_ELEVATION_ABOVE_FNAM. 
"""
function summit_markers(sb::SheetBuilder)
    promlev_prom = get_config_value("Markers", "Prominence level [m], prominent summit", Int; nothing_if_not_found = false)
    promlev_obsc = get_config_value("Markers", "Prominence level [m], obscure summit", Int; nothing_if_not_found = false)
    symbol_prom = get_config_value("Markers", "Symbol prominent summit", String; nothing_if_not_found = false)
    symbol_obsc = get_config_value("Markers", "Symbol obscure summit", String; nothing_if_not_found = false)
    symbol_size_prom = get_config_value("Markers", "Size prominent summit symbol", Int; nothing_if_not_found = false)
    symbol_size_obsc = get_config_value("Markers", "Size obscure summit symbol", Int; nothing_if_not_found = false)
    prom_levels = [promlev_obsc, promlev_prom]
    symbols = [symbol_obsc, symbol_prom]
    symbol_sizes = [symbol_size_obsc, symbol_size_prom]
    summit_markers(full_folder_path(sb), sb.cell_iter, cell_to_utm_factor(sb), sb.f_I_to_utm, prom_levels, symbols, symbol_sizes)
end
function summit_markers(fofo, cell_iter, cell2utm, f_I_to_utm, prom_levels, symbols, symbol_sizes)
    if isfile(joinpath(fofo, COMPOSITE_FNAM))
        @debug "    $COMPOSITE_FNAM in $fofo already exists. Exiting `summit_markers`"
        return true
    end
    if isfile(joinpath(fofo, MARKERS_FNAM))
        @debug "    $MARKERS_FNAM in $fofo already exists. Exiting `summit_markers`"
        return true
    end
    if ! isfile(joinpath(fofo, CONSOLIDATED_FNAM))
        @debug "    $CONSOLIDATED_FNAM in $fofo does not exist. Exiting `ridge_overlay`"
        return false
    end
    bwres = _summit_markers(fofo, cell_iter, cell2utm, f_I_to_utm, prom_levels, symbols, symbol_sizes)
    # Feedback
    display_if_vscode(bwres)
    # Save
    ffna = joinpath(fofo, MARKERS_FNAM)
    @debug "    Saving $ffna"
    save_png_with_phys(ffna, map(bwres) do pix
        pix == true && return RGBA{N0f8}(0, 0, 0, 1)
        RGBA{N0f8}(0., 0, 0, 0)
    end)
    true
end

"""
    _summit_markers(fofo, cell_iter, cell2utm, f_I_to_utm, prom_levels, symbols, symbol_sizes)
    --> Image

- Identify 'prominent' and 'obscure' summits.
- Summarize in .csv file for potential name tagging: SUMMITS_FNAM.
- Also writes the intermediate matrix 'maximum elevation above' to MAX_ELEVATION_ABOVE_FNAM. 
- Mark symbols in output image. 
"""
function _summit_markers(fofo, cell_iter, cell2utm, f_I_to_utm, prom_levels, symbols, symbol_sizes)
    ny, nx = size(cell_iter)
    si = CartesianIndices((1:cell2utm:(nx  * cell2utm), 1:cell2utm:(ny * cell2utm)))
    g = readclose(joinpath(fofo, CONSOLIDATED_FNAM))
    z = transpose(g.A[si])
    prom = find_prominence_and_write_max_elevation_above(z, joinpath(fofo, MAX_ELEVATION_ABOVE_FNAM))
    write_prominence_to_csv(prom, round.(z), f_I_to_utm, cell_iter, prom_levels, joinpath(fofo, SUMMITS_FNAM))
    @debug "    Draw summit marks"
    draw_summit_marks(prom, prom_levels, symbols, symbol_sizes)
end


"""
    find_prominence_and_write_max_elevation_above(z, ffna)
    --> Matrix

The output contains prominence values pertinent to z. Where prominence is irrelevant, NaN or zero values.
A summit is condensed to a single value, although this may fail in some cases. 

Also writes the intermediate matrix 'maximum elevation above' to MAX_ELEVATION_ABOVE_FNAM. 
"""
function find_prominence_and_write_max_elevation_above(z, ffna)
    @debug "    Find max tree "
    # We reduce the elevation resolution to whole meters in order to avoid oversegmentation,
    # which might lead to many uninteresting 'twin peaks' and heavier calculations. Many peaks are actually
    # split anyway.
    maxtree = MaxTree(round.(z))
    if isfile(ffna)
        @debug "    Load 'max elevation above' from file"
        mea = readdlm(ffna)
    else
        @debug "    Find and save 'max elevation above'"
        mea = maximum_elevation_above(z; maxtree)
        # Save
        writedlm(ffna, mea)
        open(ffna, "w") do io
            writedlm(io, mea)
        end
    end
    @debug "    Identify precise summit indices "
    summit_indices = unique_summit_indices(z, maxtree)
    @debug "    Find prominence"
    prominence(round.(z); maxtree, mea, summit_indices)
end

"""
    function prominence(z; 
        maxtree = MaxTree(z), mea = maximum_elevation_above(z; maxtree),
        summit_indices = leaf_indices(maxtree))
    ---> Matrix

The output contains prominence values pertinent to z, for
the indices supplied by `summit_indices`.
"""
function prominence(z; 
    maxtree = MaxTree(z), mea = maximum_elevation_above(z; maxtree),
    summit_indices = leaf_indices(maxtree))
    #
    # Check arguments
    maxtree.rev && throw(ArgumentError("maxtree was created with 'reverse' option."))
    # Allocate result matrix
    promin = similar(z)
    fill!(promin, NaN)
    for i in summit_indices
        summit_elevation = z[i]
        # Recursively walk down from the leaf node until meeting a dominant (taller) 
        # summit's zone. Return the meeting index.
        i_lp = lowest_parent_in_summit_zone(i, mea, maxtree, summit_elevation)
        this_prominence = summit_elevation - z[i_lp]
        promin[i] = this_prominence
    end
    promin
end

function write_prominence_to_csv(prom, z, f_I_to_utm, cell_iter, prom_levels, filename)
    vI = filter(cell_iter) do I
        p = prom[I]
        ! isnan(p) && p >= minimum(prom_levels)
    end
    vpr = prom[vI]
    vz = z[vI]
    # Utm positions
    vutm = map(f_I_to_utm, vI)
    # Sheet indices
    vi = map(I -> I.I, vI)
    # Order by prominence
    order = sortperm(vpr; rev = true)
    # Nest and prepare for output
    vectors = [vpr[order], vz[order], vutm[order], vi[order]]
    headers = ["Prominence_m", "Elevation_m", "Utm", "Sheet_index"]
    widths = [20, 20, 20, 20]  
    write_vectors_to_csv(filename, headers, vectors, widths)
end


function draw_summit_marks(prominence, prom_levels, symbols, symbol_sizes)
    length(prom_levels) == length(symbols) == length(symbol_sizes) == 2 || throw(ArgumentError("Length not 2"))
    # The prominence table is assumed oriented like an image. Equilateral triangles
    # are drawn with a horizontal line at bottom.
    #
    # Allocate black-and white image.
    img = zeros(Gray{Bool}, size(prominence)...)
    # Unpack arguments
    p1, p2 = prom_levels
    p1 < p2 || throw(ArgumentError("prom_level $prom_levels"))
    sy1, sy2 = symbols
    si1, si2 = symbol_sizes
    #
    prominent_indices = filter(i -> prominence[i] >= p2, CartesianIndices(prominence))
    obscure_indices = filter( i -> p1 <= prominence[i] < p2, CartesianIndices(prominence))
    # <-- TEMP: Since we do not yet cross-check elevations above with the neighbouring sheets,
    # we do not trust that indicated summits on the sheet border are actual summits.
    R_internal = CartesianIndices((2:size(prominence)[1] - 1, 2:size(prominence)[2] - 1))
    prominent_indices_nob = filter(I -> I ∈ R_internal, prominent_indices)
    obscure_indices_nob = filter(I -> I ∈ R_internal, obscure_indices)
    # TEMP -->
    mark_at!(img, obscure_indices_nob, si1, sy1)
    mark_at!(img, prominent_indices_nob, si2, sy2)
end


function maximum_elevation_above(z; maxtree = MaxTree(z))
    # Allocate result matrix
    mea = similar(z)
    fill!(mea, NaN)
    # We're starting with leaf nodes (i.e. the maxima)
    # and recursively visit all parents of it with 
    # the information about the tallness of it's tallest
    # descendants. The recursion returns when
    #    - a  node is reached
    #    - or we reach a node that already knows about a taller parent.
    # At return, we pick the next node in this loop.
    vi = leaf_indices(maxtree)
    for i in vi
        # The leaf is most likely the parent of a 'reference node'.
        # That 'reference node' is most likely NOT a summit.
        leaf_elevation = round(z[i])
        parent_i = maxtree.parentindices[i]
        # Call the recursive, in_place function.
        _mea!(mea, maxtree, leaf_elevation, parent_i)
    end
    # At this point, all leaf nodes still contain just NaN.
    # Now copy the value from each leaf's reference node.
    for i in vi
        parent_i = maxtree.parentindices[i]
        mea[i] = mea[parent_i]
    end
    mea
end

function _mea!(mea, maxtree, leaf_elevation, i)
    previous_maximum_elevation_above = mea[i]
    if isnan(previous_maximum_elevation_above)
        # Remember this leaf. Most likely, another and taller
        # leaf will overwrite this later.
        mea[i] = leaf_elevation
    elseif previous_maximum_elevation_above >= leaf_elevation
        # Just leave quietly. Stop the recursion here.
        return mea
    else 
        # Forget about what we thought was the summit elevation.
        # Remember this instead!
        mea[i] = leaf_elevation
    end
    # Let's look for our parent.
    parent_i = maxtree.parentindices[i]
    if parent_i == i
        # ⬆ We are our own parent, the root, lowest level in the map.
        # Stop the recursion here.
        return mea
    end
    # ⬇ Recurse: Tell current parent (which is lower down) 
    # about the leaf's elevation. It's the actual summit above us for all we know.
    _mea!(mea, maxtree, leaf_elevation, parent_i)
end

function lowest_parent_in_summit_zone(i, mea, maxtree, summit_elevation)
    ic = CartesianIndices(mea)[i].I
    max_elevation_above = mea[i]
    if isnan(max_elevation_above)
        throw(ArgumentError("The maximum_elevation_above matrix has a NaN value at $(ic)"))
    elseif max_elevation_above > summit_elevation
        # ⬆ Reached down to a taller elevation's summit zone
        # Stop the recursion here.
        return i
    end
    parent_i = maxtree.parentindices[i]
    if parent_i == i
        # ⬆ We are our own parent, the root, lowest level in the map.
        # Stop the recursion here.
        return i
    end
    # ⬇ still in the same summit zone. Proceed downslope.
    # Recurse
    lowest_parent_in_summit_zone(parent_i, mea, maxtree, summit_elevation)
end

function leaf_indices(maxtree; parent_set = Set(maxtree.parentindices))
    # maxtree.parentindices is a matrix with mostly 
    # repeating parents. That means, a minority of nodes are parents.
    # Most nodes are not referred, and are thus considered leaves.
    unique(filter(i -> i ∉ parent_set, eachindex(maxtree.parentindices)))
end

function parent_of_leaf_indices(maxtree; parent_set = Set(maxtree.parentindices))
    li = leaf_indices(maxtree; parent_set)
    unique(maxtree.parentindices[li])
end


"""
    core_family_dictionary(maxtree)
    --> Dict{Int64, Vector{Int64}}

Dictionary of components where a parent (key) 
  - has children (value) but not grandchildren.
  - is not a root.
"""
function core_family_dictionary(maxtree)
    # All the parent indices, but just one of each. 
    parent_set = Set(maxtree.parentindices)
    spli = Set(parent_of_leaf_indices(maxtree; parent_set))
    # Build the family dictionary: parent => [child1, child2..]
    dic = Dict{Int, Vector{Int}}()
    for (i, parent_i) in enumerate(maxtree.parentindices)
        if parent_i ∈ spli && parent_i !== i
            dic[parent_i] = push!(get(dic, parent_i, Int[]), i)
        end
    end
    # Throw out the entire families where a child is also a parent.
    filter(dic) do (_, children)
        isempty(intersect(parent_set, children))
    end
end


"""
    unique_summit_indices(z, maxtree)
    -- Set{Int64}


The maxtree has been derived from matrix z (a rounded version of it).

We can imagine several types of summits as represented by the maxtree.

Case A: 
    A local maximum has no neigbouring elevations that are identical. This
    will be a leaf node, e.g. '2' here:

    1 1 1
    1 2 1
    1 1 1

    Case A returns Set(5), cartesian index [2, 2]


 Case B: 
    A local maximum has neigbouring elevations that are identical, 
    or at least identical when rounded. Some rounding avoids
    oversegmentation and simplifies calculations. Example of values
    before rounding which might lead to this:

      1.0  1.0  1.0
      1.0  2.0  2.2
      1.0  2.1  1.0

   After rounding and taking the maxtree, the '2.2' and '2.1' would 
   become leaf nodes, and 2.0 a reference node.

   `unique_summit_indices(B, MaxTree(round.(B))` 
   --> Set(8) , i.e. cartesian [2, 3]



Another example of case B returns Set(7), cartesian [3, 2]:

    21.0  19.0  19.0  16.0
    22.0  22.0  23.0  22.0
    22.0  23.0  23.0  21.0
    22.0  22.0  22.0  16.0

"""
function unique_summit_indices(z, maxtree)
    # Note:  The number of candidates for summits is typically 
    # in the millions for one sheet. Hence, we avoid nested loops.
    #
    # Check arguments
    maxtree.rev && throw(ArgumentError("maxtree was created with 'reverse' option."))
    # Dictionary of components where a parent has children but not grandchildren
    dic_fam = core_family_dictionary(maxtree)
    summits = Set{Int}()
    # For each family, add the index of it's tallest member, 
    # be it a parent (reference node) or child (leaf node).
    for (parent, children) in dic_fam
        z_parent = z[parent]
        tallest_z, tallest_child_no = findmax(z[children])
        if tallest_z > z_parent
            # A child is tallest in this family.
            push!(summits, children[tallest_child_no])
        else
            # The parent is tallest in this family.
            push!(summits, parent)
        end
    end
    summits
end
