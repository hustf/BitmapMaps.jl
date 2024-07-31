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
    summit_prominence_sheet_internal(sb::SheetBuilder)
"""
function summit_prominence_sheet_internal(sb::SheetBuilder)
    promlev_prom = get_config_value("Markers", "Prominence level [m], prominent summit", Int; nothing_if_not_found = false)
    promlev_obsc = get_config_value("Markers", "Prominence level [m], obscure summit", Int; nothing_if_not_found = false)
    symbol_prom = get_config_value("Markers", "Symbol prominent summit", String; nothing_if_not_found = false)
    symbol_obsc = get_config_value("Markers", "Symbol obscure summit", String; nothing_if_not_found = false)
    symbol_size_prom = get_config_value("Markers", "Size prominent summit symbol", Int; nothing_if_not_found = false)
    symbol_size_obsc = get_config_value("Markers", "Size obscure summit symbol", Int; nothing_if_not_found = false)
    prom_levels = [promlev_obsc, promlev_prom]
    symbols = [symbol_obsc, symbol_prom]
    symbol_sizes = [symbol_size_obsc, symbol_size_prom]
    summit_prominence_sheet_internal(full_folder_path(sb), sb.cell_iter, cell_to_utm_factor(sb), sb.f_I_to_utm, prom_levels, symbols, symbol_sizes)
end

function summit_prominence_sheet_internal(fofo, cell_iter, cell2utm, f_I_to_utm, prom_levels, symbols, symbol_sizes)
    if isfile(joinpath(fofo, COMPOSITE_FNAM))
        @debug "    $COMPOSITE_FNAM in $fofo already exists. Exiting `summit_prominence_sheet_internal`"
        return true
    end
    if isfile(joinpath(fofo, MARKERS_FNAM))
        @debug "    $MARKERS_FNAM in $fofo already exists. Exiting `summit_prominence_sheet_internal`"
        return true
    end
    bwimg = _summit_prominence_sheet_internal(fofo, cell_iter, cell2utm, f_I_to_utm, prom_levels, symbols, symbol_sizes)
    # Feedback
    display_if_vscode(bwimg)
    # Save
    ffna = joinpath(fofo, MARKERS_FNAM)
    @debug "    Saving $ffna"
    save_png_with_phys(ffna, map(bwimg) do pix
        pix == true && return RGBA{N0f8}(0, 0, 0, 1)
        RGBA{N0f8}(0., 0, 0, 0)
    end)
    true
end

function _summit_prominence_sheet_internal(fofo, cell_iter, cell2utm, f_I_to_utm, prom_levels, symbols, symbol_sizes)
    # Get elevation matrix. This samples every point regardless of cell_to_utm_factor
    g = readclose(joinpath(fofo, CONSOLIDATED_FNAM))
    ny, nx = size(cell_iter)
    source_indices = (1:cell2utm:(nx  * cell2utm), 1:cell2utm:(ny * cell2utm))
    si = CartesianIndices(source_indices)
    # For largish sheets, we need to cut down on the number of 
    # levels, otherwise we may get stack overflow errors.
    # Look for stack traces with e.g. _mea!...(repeats 52180 times)
    prom = find_prominence_and_write_max_elevation_above(Float16.(g.A[si]), fofo)
    write_prominence_location(prom, f_I_to_utm)
    @debug "    Draw summit marks"
    draw_summit_marks(transpose(prom), prom_levels, symbols, symbol_sizes)
end

function find_prominence_and_write_max_elevation_above(z, fofo)
    @debug "    Find max tree "
    maxtree = MaxTree(z)
    @debug "    Find max elevation above"
    mea = maximum_elevation_above(z; maxtree)
    # Save
    ffna = joinpath(fofo, MAX_ELEVATION_ABOVE_FNAM)
    @debug "    Saving $ffna"
    #save_png_with_phys(ffna, Gray.(mea))
    @debug "    Find prominence"
    prominence(z; maxtree, mea )
end

function write_prominence_location(prom, f_I_to_utm)
    @info "write_prominence_location not implemented. Consider how much time it would take to 
           recalculate given the file maximum elevation above."
end

function draw_summit_marks(prominence, prom_levels, symbols, symbol_sizes)
    length(prom_levels) == length(symbols) == length(symbol_sizes) == 2 || throw(ArgumentError("Length not 2"))

    # The prominence table is assumed oriented like an image. Equilateral triangles
    # are drawn with a horizontal line at bottom.
    #
    # Allocate black-and white image.
    img = zeros(Gray{Bool}, size(prominence)...)
    #
    p1, p2 = prom_levels
    p1 < p2 || throw(ArgumentError("prom_level $prom_levels"))
    sy1, sy2 = symbols
    si1, si2 = symbol_sizes
    #
    prominent_indices = filter(i -> prominence[i] >= p2, CartesianIndices(prominence))
    obscure_indices = filter( i -> p1 <= prominence[i] < p2, CartesianIndices(prominence))
    mark_at!(img, obscure_indices, si1, sy1)
    mark_at!(img, prominent_indices, si2, sy2)
end


function prominence(z; 
    maxtree = MaxTree(z), mea = maximum_elevation_above(z; maxtree))
    # Allocate result matrix
    promin = similar(z)
    fill!(promin, NaN)
    # Identify (part of) all summits as leaf nodes.
    # TODO bugfix
    vi = leaves(maxtree)
    for (count, i) in enumerate(vi)
        summit_elevation = z[i]
        ic = CartesianIndices(z)[i].I
        #@debug "$count of $(length(vi)) → $ic  Summit elevation $(Float64(summit_elevation))"
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
        #@warn "Finished summit leaf $ic, found it has prominence $(promin[i]) "
    end
    promin
end



function maximum_elevation_above(z; maxtree = MaxTree(z))
    # Allocate result matrix
    mea = similar(z)
    fill!(mea, NaN)
    # We're starting with leaf nodes (i.e. the maxima,
    # or parents of none)
    # and recursively visit all parents of it with 
    # the information about the tallness of it's tallest
    # descendants. The recursion returns when
    #    - a  node is reached
    #    - or we reach a node that already knows about a taller parent.
    # At return, we pick the next node in this loop.
    vi = leaves(maxtree)
    for (count, i) in enumerate(vi)
        summit_elevation = z[i]
        # Call the recursive, in_place function.
        ic = CartesianIndices(mea)[i].I
        #@debug "$count of $(length(vi)) → $ic  Summit elevation $(Float64(summit_elevation))"
        _mea!(mea, maxtree, summit_elevation, i)
    end
    mea
end
function _mea!(mea, maxtree, summit_elevation, i)
    ic = CartesianIndices(mea)[i].I
    previous_maximum_elevation_above = mea[i]
    if isnan(previous_maximum_elevation_above)
        mea[i] = summit_elevation
    elseif previous_maximum_elevation_above >= summit_elevation
        #@debug "⬆ $ic  Previous >= $(Float64(summit_elevation))"
        return mea
    else 
        mea[i] = summit_elevation
    end
    parent_i = maxtree.parentindices[i]
    if parent_i == i
        #@debug "⬆ $ic Is the root."
        return mea
    end
    # Recurse
#    @debug "⬇ $ic  Tell parent $(CartesianIndices(mea)[parent_i].I)"
    _mea!(mea, maxtree, summit_elevation, parent_i)
end

function lowest_parent_in_summit_zone(i, mea, maxtree, summit_elevation)
    ic = CartesianIndices(mea)[i].I
    max_elevation_above = mea[i]
    if isnan(max_elevation_above)
        throw(ArgumentError("The maximum_elevation_above matrix has a NaN value at $(ic)"))
    elseif max_elevation_above > summit_elevation
        #@debug """⬆ $ic  Reached a taller elevation's summit zone at $ic:
        #    $(max_elevation_above) > $(Float64(summit_elevation))"""
        return i
    end
    parent_i = maxtree.parentindices[i]
    if parent_i == i
        #@debug "⬆ $ic Is the root."
        return i
    end
    #@debug "⬇ $ic  is in the same summit zone $(Float64(summit_elevation)). Proceed downslope."
    # Recurse
    lowest_parent_in_summit_zone(parent_i, mea, maxtree, summit_elevation)
end

#haschild(i, maxtree) = i ∈ maxtree.parentindices
#leaves(maxtree) = filter(i-> ! haschild(i, maxtree), eachindex(maxtree.parentindices))
function leaves(maxtree)
    parent_set = Set(maxtree.parentindices)
    return filter(i -> i ∉ parent_set, eachindex(maxtree.parentindices))
end