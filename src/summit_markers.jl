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

function _summit_markers(fofo, cell_iter, cell2utm, f_I_to_utm, prom_levels, symbols, symbol_sizes)
    ny, nx = size(cell_iter)
    si = CartesianIndices((1:cell2utm:(nx  * cell2utm), 1:cell2utm:(ny * cell2utm)))
    g = readclose(joinpath(fofo, CONSOLIDATED_FNAM))

    z = transpose(g.A[si])
    prom = find_prominence_and_write_max_elevation_above(z, joinpath(fofo, MAX_ELEVATION_ABOVE_FNAM))
    write_prominence_data_to_csv(prom, round.(z), f_I_to_utm, cell_iter, prom_levels, joinpath(fofo, SUMMITS_FNAM))
    @debug "    Draw summit marks"
    draw_summit_marks(prom, prom_levels, symbols, symbol_sizes)
end

function find_prominence_and_write_max_elevation_above(z, ffna)
    @debug "    Find max tree "
    # We reduce the elevation resolution to whole meters in order to avoid oversegmentation,
    # which migt lead to many uninteresting 'twin peaks' and heavier calculations. Many peaks are actually
    # split anyway.
    maxtree = MaxTree(round.(z))
    if isfile(ffna)
        @debug "    Find max elevation above from file"
        mea = readdlm(ffna)
    else
        @debug "    Find max elevation above"
        mea = maximum_elevation_above(round.(z); maxtree)
            # Save
        @debug "    Saving $ffna"
        writedlm(ffna, mea)
        open(ffna, "w") do io
            writedlm(io, mea)
        end
    end
    @debug "    Find prominence"
    prominence_clusters(round.(z); maxtree, mea)
    # TODO: Remove clusters, pick in each cluster from the unrounded elevation.
end

function write_prominence_data_to_csv(prom, z, f_I_to_utm, cell_iter, prom_levels, filename)
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
   
function prominence_clusters(z; 
    maxtree = MaxTree(z), mea = maximum_elevation_above(z; maxtree))
    # Check argument
    maxtree.rev && throw(ArgumentError("maxtree was created with 'reverse' option."))
    # Allocate result matrix
    promin = similar(z)
    fill!(promin, NaN)
    # Identify (part of) all summits.
    vi = leaf_indices(maxtree)

    for i in vi
        summit_elevation = z[i]
        # Recursively walk down from the leaf node until meeting a dominant (taller) 
        # summit's zone. Return the meeting index.
        i_lp = lowest_parent_in_summit_zone(i, mea, maxtree, summit_elevation)
        this_prominence = summit_elevation - z[i_lp]
        # Look for potential reference nodes:
        # All summits include leaf nodes, but if 
        # the summit is flat, one reference node is a parent 
        # of leaf nodes on the summit. 
        i_parent = maxtree.parentindices[i]
        if z[i_parent] == summit_elevation
            # Propagate leaf prominence to the reference node
            promin[i_parent] = this_prominence
        end
        promin[i] = this_prominence
        # Finished summit /leaf, found its prominence
    end
    promin
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
        leaf_elevation = z[i]
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

function leaf_indices(maxtree)
    # maxtree.parentindices is a matrix with mostly 
    # repeating parents. That is, just a few nodes are parents.
    parent_set = Set(maxtree.parentindices)
    # Most nodes are not referred, and are thus considered leaves.
    return filter(i -> i ∉ parent_set, eachindex(maxtree.parentindices))
end

