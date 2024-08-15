# This file concerns finding summit prominence.
# The prominence of some summits actually depends
# on geography which lay outside the limits of one sheet.
# We can't fit every sheet in memory at once, so 
# a sheet contact solution should be found after al the initial single sheet
# calculations are finished. 
# This file prepares for the sheet contact by:
# - prepare and save 'maximum elevation above' (mea) table.
# -  export the outermost row or column, depending on direction,
#      to the neighouring sheets. We export elevation (z) and 'maximum elevation' as
#      separate files.


"""
    summit_markers(sb::SheetBuilder)
    summit_markers(fofo, cell_iter, cell2utm, f_I_to_utm, prom_levels, symbols, symbol_sizes, dic_neighbour)
    --> Bool

- Identify 'prominent' and 'obscure' summits.
- Mark symbols in image MARKERS_FNAM. 
- Summarize in .csv file for potential name tagging: SUMMITS_FNAM.
- Also writes the intermediate matrix 'maximum elevation above' to MAX_ELEVATION_ABOVE_FNAM. 
"""
function summit_markers(sb::SheetBuilder)
    # Harvest info
    promlev_prom = get_config_value("Markers", "Prominence level [m], prominent summit", Int; nothing_if_not_found = false)
    promlev_obsc = get_config_value("Markers", "Prominence level [m], obscure summit", Int; nothing_if_not_found = false)
    symbol_prom = get_config_value("Markers", "Symbol prominent summit", String; nothing_if_not_found = false)
    symbol_obsc = get_config_value("Markers", "Symbol obscure summit", String; nothing_if_not_found = false)
    symbol_size_prom = get_config_value("Markers", "Size prominent summit symbol", Int; nothing_if_not_found = false)
    symbol_size_obsc = get_config_value("Markers", "Size obscure summit symbol", Int; nothing_if_not_found = false)
    prom_levels = [promlev_obsc, promlev_prom]
    symbols = [symbol_obsc, symbol_prom]
    symbol_sizes = [symbol_size_obsc, symbol_size_prom]
    # The neighbouring folders are perhaps being created at this point.
    # Postponing collecting the paths to neighbours might save calculation iterations, but also increase complexity.
    # So we do it early.
    dic_neighbour = neighbour_folder_dict(sb)
    summit_markers(full_folder_path(sb), sb.cell_iter, cell_to_utm_factor(sb), sb.f_I_to_utm, prom_levels, symbols, symbol_sizes, dic_neighbour)
end
function summit_markers(fofo, cell_iter, cell2utm, f_I_to_utm, prom_levels, symbols, symbol_sizes, dic_neighbour)
    # Early exits: 
    if ! isfile(joinpath(fofo, CONSOLIDATED_FNAM))
        @debug "    $CONSOLIDATED_FNAM in $fofo does not exist. Exiting `summit_markers`"
        return false
    end

    #= TEMP  TODO - we should have some criterion when more iteration is not necessary.
    if isfile(joinpath(fofo, MARKERS_FNAM))
        @debug "    $MARKERS_FNAM in $fofo already exists. Exiting `summit_markers`"
        return true
    end
    =#
    if isfile(joinpath(fofo, COMPOSITE_FNAM))
        @debug "    $COMPOSITE_FNAM in $fofo is being deleted, because we iterate for sheet iteration on `summit_markers`."
        rm(joinpath(fofo, COMPOSITE_FNAM))
    end
    bwres = _summit_markers(fofo, cell_iter, cell2utm, f_I_to_utm, prom_levels, symbols, symbol_sizes, dic_neighbour)
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
    _summit_markers(fofo, cell_iter, cell2utm, f_I_to_utm, prom_levels, symbols, symbol_sizes, dic_neighbour)
    --> Image

- Identify 'prominent' and 'obscure' summits.
- Summarize in .csv file for potential name tagging: SUMMITS_FNAM.
- Also writes the intermediate matrix 'maximum elevation above' to MAX_ELEVATION_ABOVE_FNAM. 
- Mark symbols in output image. 
"""
function _summit_markers(fofo, cell_iter, cell2utm, f_I_to_utm, prom_levels, symbols, symbol_sizes, dic_neighbour)
    ny, nx = size(cell_iter)
    si = CartesianIndices((1:cell2utm:(nx  * cell2utm), 1:cell2utm:(ny * cell2utm)))
    g = readclose(joinpath(fofo, CONSOLIDATED_FNAM))
    z = transpose(g.A[si])
    prom = find_prominence_and_write_max_elevation_above(z, fofo; dic_neighbour)
    # Prom is a large matrix with mostly NaN values. 
    # Let's find the relevant indices.
    # Not interesting: 
    #   1) is NaN
    #   2) prominence < defined in the .ini file
    #   3) elevation < 200, hardcoded. This removes a lot of power-line "summits", waves and other artifacts.
    vI = filter(cell_iter) do I
        p = prom[I]
        ! isnan(p) && p >= minimum(prom_levels) && z[I] > 200
    end
    write_prominence_to_csv(prom, round.(z), f_I_to_utm, vI, joinpath(fofo, SUMMITS_FNAM))
    draw_summit_marks(prom, vI, prom_levels, symbols, symbol_sizes)
end


"""
    find_prominence_and_write_max_elevation_above(z, fo; dic_neighbour = Dict{Symbol, String}())
    --> Matrix

Writes the result to MAX_ELEVATION_ABOVE_FNAM.

Also distributes boundary conditions, eight files with a vector each, to neighbouring sheet's folders 
as specified in `dic_neighbour`.

The output contains prominence values pertinent to z. Where prominence is irrelevant, NaN or zero values.
A summit is condensed to a single value, although this may fail in some cases. 

Also writes the intermediate matrix 'maximum elevation above' to MAX_ELEVATION_ABOVE_FNAM. 
"""
function find_prominence_and_write_max_elevation_above(z, fo; dic_neighbour = Dict{Symbol, String}())
    @assert isdir(fo)
    @debug "    Find max tree "
    # We reduce the elevation resolution to whole meters in order to avoid oversegmentation,
    # which might lead to many uninteresting 'twin peaks' and heavier calculations. Many peaks are actually
    # split anyway.
    maxtree = MaxTree(round.(z))
    # Find the most important indices
    summit_indices = distinct_summit_indices(z, maxtree)
    # Retrieve boundary conditions
    boundaries = read_boundaries_z_mea(fo, size(z, 1), size(z, 2))
    #
    ffna = joinpath(fo, MAX_ELEVATION_ABOVE_FNAM)
    # Recalculate or load?
    #if isfile(ffna)
    #    @debug "    Load 'max elevation above' from file"
        # TEMP consider recalculation because of changed boundary conditions
    #    mea = eltype(z).(readdlm(ffna))
    #else
        mea = maximum_elevation_above(z; maxtree, summit_indices, boundaries)
        # Save
        open(ffna, "w") do io
            writedlm(io, mea)
        end
    #end
    @debug "    Distribute boundary conditions"
    distribute_boundary_conditions(dic_neighbour, z, mea)
    @debug "    Find prominence"
    prominence(z, summit_indices, mea; maxtree)
end

"""
    prominence(z, summit_indices, mea; maxtree = MaxTree(z))
    ---> Matrix

Outputs [prominence](https://en.wikipedia.org/wiki/Topographic_prominence) values pertinent to z, for
the indices supplied by `summit_indices`.

# Example

```
julia> z = [ 0.0   1.0   0.0  -3.0   -8.0                                                                                                                                              
            1.0   2.0   1.0  -2.0   -7.0                                                                                                                                               
            0.0   1.0   0.0   4.0    3.0                                                                                                                                               
           -3.0  -2.0  -3.0  -6.0  -11.0];

julia> maxtree = MaxTree(round.(z));

julia> summit_indices = distinct_summit_indices(z, maxtree);

julia> mea = maximum_elevation_above(z; maxtree, summit_indices);

julia> prominence(z, summit_indices, mea; maxtree)
4×5 Matrix{Float64}:
 NaN  NaN    NaN  NaN    NaN
 NaN    2.0  NaN  NaN    NaN
 NaN  NaN    NaN   15.0  NaN
 NaN  NaN    NaN  NaN    NaN

```
"""
function prominence(z, summit_indices, mea; maxtree = MaxTree(z))
    # Check arguments
    maxtree.rev && throw(ArgumentError("maxtree was created with 'reverse' option."))
    # Define valid regions for prominence values
    R = CartesianIndices(z)
    R_internal = CartesianIndices((2 : size(z, 1) - 1, 2 : size(z, 2) - 1))
    # Allocate result matrix
    promin = similar(z)
    fill!(promin, NaN)
    for i in R[collect(summit_indices)]
        if i ∈ R_internal
            # Even though we consider boundary conditions,
            # it requires some work to define those for the external boundaries
            # to a collection of sheets. At the same time, a summit on the edges of 
            # a sheet is extremely unlikely. Everything considered, we discard all such summits.
            summit_elevation = z[i]
            # Recursively walk down from the leaf node until meeting a dominant (taller) 
            # summit's zone. Return the meeting index.
            i_lp = lowest_parent_in_summit_zone(i, mea, maxtree, summit_elevation)
            this_prominence = summit_elevation - z[i_lp]
            promin[i] = this_prominence
        end
    end
    promin
end

"""
    write_prominence_to_csv(prom, z, f_I_to_utm, cell_iter, filename)
    ---> Nothing

Writes elevation, prominence, position and index.
"""
function write_prominence_to_csv(prom, z, f_I_to_utm, indices, filename)
    @debug "    Saving $filename"
    # Interesting values
    vpr = prom[indices]
    vz = z[indices]
    # Utm positions
    vutm = map(f_I_to_utm, indices)
    # Sheet indices
    vi = map(I -> I.I, indices)
    # Order by elevation
    order = sortperm(vz; rev = true)
    # Nest and prepare for output
    vectors = [vz[order], vpr[order], vutm[order], vi[order]]
    headers = ["Elevation_m", "Prominence_m", "Utm", "Sheet_index"]
    widths = [20, 20, 20, 20]  
    write_vectors_to_csv(filename, headers, vectors, widths)
end


function draw_summit_marks(prominence, indices, prom_levels, symbols, symbol_sizes)
    length(prom_levels) == length(symbols) == length(symbol_sizes) == 2 || throw(ArgumentError("Length not 2"))
    #
    @debug "    Draw summit marks"
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
    prominent_indices = filter(indices) do i
        prominence[i] >= p2
    end
    obscure_indices = filter(indices) do i
        p1 <= prominence[i] < p2 >= p2
    end
    mark_at!(img, obscure_indices, si1, sy1)
    mark_at!(img, prominent_indices, si2, sy2)
end

"""
    maximum_elevation_above(z::T; 
        maxtree = MaxTree(z), 
        summit_indices = distinct_summit_indices(z, maxtree)) where T <: AbstractArray
    ---> T

Finding the maxtree from unrounded values of z (as is the default with no keywords supplied)
is slow. A weak recommendation is to supply both keywords, with `maxtree = MaxTree(round.(z))`.

See 'prominence'.
"""
function maximum_elevation_above(z::T; 
    maxtree = MaxTree(z), 
    summit_indices = distinct_summit_indices(z, maxtree),
    boundaries = nothing) where T <: AbstractArray
    #

    # Check arguments
    eltype(z) <: AbstractFloat || throw(ArgumentError("eltype(z) = $(eltype(z)) is not <: AbstractFloat"))
    # This algorithm should be more effective if we start with the tallest summits.
    # Consider that summit_indices is a Set. Make a sorted vector.
    vsummit_indices = let v = collect(summit_indices)
        v[sortperm(z[v]; rev = true)]
    end
    # Allocate result matrix
    mea = fill(eltype(z)(NaN), size(z))
    # Visit direct and indicert parents of summits with 
    # the information about the tallness of its tallest
    # descendants. The recursion returns when
    #    - a root is reached
    #    - or we reach a node that already knows about a taller parent.
    # At return, we pick the next node in this loop.
    if isnothing(boundaries) || all(sum.(boundaries) .== 0)
        # This will be run if none of the neighouring sheets
        # have calculated their 'mea' matrices yet.
        @debug "    Find 'max elevation above' for this sheet in isolation"
        for i in vsummit_indices
            max_elevation = z[i]
            parent_i = maxtree.parentindices[i]
            # Call the recursive, in_place function.
            _mea!(mea, maxtree, max_elevation, parent_i)
        end
    else
        # This takes 'mea' from neighbouring sheets into account.
        @debug "    Find 'max elevation above', including boundary conditions from other sheets."
        f! = func_mea_contact!(maxtree, z, boundaries)
        for i in vsummit_indices
            max_elevation = z[i] # TEMP. If the summit lies on a boundary, this is incorrect...
            parent_i = maxtree.parentindices[i]
            # Call the recursive, in_place function.
            f!(mea, max_elevation, parent_i)
        end
    end
    # At this point, summit indices and some of their potential leaf nodes still contain just NaN.
    # Now copy the value from their parents.
    for i in vsummit_indices
        parent_i = maxtree.parentindices[i]
        mea[i] = mea[parent_i]
    end
    for i in leaf_indices(maxtree)
        parent_i = maxtree.parentindices[i]
        if isnan(mea[i])
            mea[i] = mea[parent_i]
        else
            if mea[i] !== mea[parent_i]
                @show mea[i] mea[parent_i]
                throw("unexpected")
            end
        end
    end
    mea
end

function _mea!(mea, maxtree, leaf_elevation, i)
    previous_max_elevation_above = mea[i]
    if isnan(previous_max_elevation_above)
        # Remember this leaf. Most likely, another and taller
        # leaf will overwrite this later.
        mea[i] = leaf_elevation
    elseif previous_max_elevation_above >= leaf_elevation
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


#=
=#
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
    distinct_summit_indices(z, maxtree)
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

   `distinct_summit_indices(B, MaxTree(round.(B))` 
   --> Set(8) , i.e. cartesian [2, 3]



Another example of case B returns Set(7), cartesian [3, 2]:

    21.0  19.0  19.0  16.0
    22.0  22.0  23.0  22.0
    22.0  23.0  23.0  21.0
    22.0  22.0  22.0  16.0

"""
function distinct_summit_indices(z, maxtree)
    # Note:  The number of candidates for summits is typically 
    # in the millions for one sheet. Hence, we avoid nested loops.
    #
    @debug "    Identify precise summit indices "
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
