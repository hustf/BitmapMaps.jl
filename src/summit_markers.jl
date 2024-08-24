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
    summit_markers(fofo, cell_iter, cell2utm, f_I_to_utm, prom_levels, symbols, symbol_sizes, min_stress, dic_neighbour)
    --> Bool

- Identify 'prominent' and 'obscure' summits.
- Mark symbols in image MARKERS_FNAM.
- Summarize in .csv file for potential name tagging: SUMMITS_FNAM.
"""
function summit_markers(sb::SheetBuilder)
    # Harvest info
    promlev_prom = get_config_value("Markers", "Prominence level [m], prominent summit", Int; nothing_if_not_found = false)
    promlev_obsc = get_config_value("Markers", "Prominence level [m], obscure summit", Int; nothing_if_not_found = false)
    symbol_prom = get_config_value("Markers", "Symbol prominent summit", String; nothing_if_not_found = false)
    symbol_obsc = get_config_value("Markers", "Symbol obscure summit", String; nothing_if_not_found = false)
    symbol_size_prom = get_config_value("Markers", "Size prominent summit symbol", Int; nothing_if_not_found = false)
    symbol_size_obsc = get_config_value("Markers", "Size obscure summit symbol", Int; nothing_if_not_found = false)
    min_stress = get_config_value("Markers", "Minimum stress level", Int; nothing_if_not_found = false)
    prom_levels = [promlev_obsc, promlev_prom]
    symbols = [symbol_obsc, symbol_prom]
    symbol_sizes = [symbol_size_obsc, symbol_size_prom]
    # The neighbouring folders are perhaps being created at this point.
    # Postponing collecting the paths to neighbours might save calculation iterations, but also increase complexity.
    # So we do it early.
    dic_neighbour = neighbour_folder_dict(sb)
    summit_markers(full_folder_path(sb), sb.cell_iter, cell_to_utm_factor(sb), sb.f_I_to_utm, prom_levels, symbols, symbol_sizes, min_stress, dic_neighbour)
end
function summit_markers(fofo, cell_iter, cell2utm, f_I_to_utm, prom_levels, symbols, symbol_sizes, min_stress, dic_neighbour)
    # Early exits:
    if ! isfile(joinpath(fofo, CONSOLIDATED_FNAM))
        @debug "    $CONSOLIDATED_FNAM in $fofo does not exist. Exiting `summit_markers`"
        return false
    end
    # Do the calculation, and saving of boundary conditions
    _summit_markers(fofo, cell_iter, cell2utm, f_I_to_utm, prom_levels, symbols, symbol_sizes, min_stress, dic_neighbour)
    true
end

"""
    _summit_markers(fofo, cell_iter, cell2utm, f_I_to_utm, prom_levels, symbols, symbol_sizes, dic_neighbour)
    --> Image

- Identify 'prominent' and 'obscure' summits.
- Summarize in .csv file: SUMMITS_FNAM.
- Mark symbols in output image.
- Modify .csv files, add name column from online lookup.
"""
function _summit_markers(fofo, cell_iter, cell2utm, f_I_to_utm, prom_levels, symbols, symbol_sizes, min_stress, dic_neighbour)
    ny, nx = size(cell_iter)
    si = CartesianIndices((1:cell2utm:(nx  * cell2utm), 1:cell2utm:(ny * cell2utm)))
    g = readclose(joinpath(fofo, CONSOLIDATED_FNAM))
    z = transpose(g.A[si])    
    prom = find_prominence(z, fofo; dic_neighbour)
    # Prom is a large matrix with mostly NaN values.
    # Let's find the relevant indices.
    # Not interesting:
    #   1) is NaN
    #   2) prominence < defined in the .ini file
    #   3) elevation < 200, hardcoded. This removes a lot of power-line "summits", waves and other artifacts.
    #   4) Artifacts > 200 (we'll do this after calculating the tensile stress equivalent.)
    minprocrit = minimum(prom_levels)
    vI = filter(cell_iter) do I
        p = prom[I]
        isnan(p) && return false
        p < minprocrit && return false
        z[I] < 200 && return false
        true
    end
    #
    # Calculate the tensile stress equvalent for summit indices
    #
    R = CartesianIndices(g.A[:,:,1])
    # Note that the Hessian components at each sample cell [0,0] is affected by region [-2:2, 2:2]
    Ω = strel(CartesianIndices((-2:2, -2:2)))
    vσ = map(vI) do I
        # Now we'd like to extract the principal tensile stress at this I. But
        # we want the stress to be consistent, disregarding downsampling.
        # So we'll take the stress level from the non-downsampled z matrix.
        # We also skip transposing it here, for assumed speed gain.
        J = CartesianIndex{2}((1 + (I.I[2] - 1) * cell2utm, (1 + (I.I[1] - 1) * cell2utm)))
        # The vincinity needed for the Hessian
        Ωⱼ = Ref(J) .+ Ω
        # If we stray outside of the entire R (unlikely), return zero as a simplification.
        any(q -> q ∉ R, Ωⱼ) && return 0.0
        #
        # The Hessian components are used to determine the principal tensile stress equivalent
        g11, _, _, g22 = hessian_components(g.A[Ωⱼ])
        # Throw away all but the central cell, return something akin to the principal tensile stress
        g11[3, 3] + g22[3,3]
    end
    #
    # We'll remove some more false peaks from this user-adjustable criterion:
    #   4) Very steep on all sides of the tallest cell, i.e. very large amplitude Hessian components.
    σII = [vσ[i] for i in 1:length(vσ) if vσ[i] >= min_stress]
    vII = [vI[i] for i in 1:length(vI) if vσ[i] >= min_stress]
    # Store the results
    waschanged = write_prominence_to_csv(prom[vII], round.(z[vII]), map(f_I_to_utm, vII), vII, σII, joinpath(fofo, SUMMITS_FNAM))
    # Add names to summits file, also save markers image
    ffna = joinpath(fofo, MARKERS_FNAM)
    if waschanged || ! isfile(ffna)
        add_names_to_csv(joinpath(fofo, SUMMITS_FNAM))
        bwres = draw_summit_marks(prom, vII, prom_levels, symbols, symbol_sizes)
        # Feedback
        display_if_vscode(bwres)
        @debug "    Saving $ffna"
        save_png_with_phys(ffna, map(bwres) do pix
            pix == true && return RGBA{N0f8}(0, 0, 0, 1)
            RGBA{N0f8}(0., 0, 0, 0)
        end)
    end
end


"""
    find_prominence(z, fofo; dic_neighbour = Dict{Symbol, String}())
    --> Matrix

Also distributes boundary conditions, eight files with a vector each, to neighbouring sheet's folders
as specified in `dic_neighbour`.

The output contains prominence values pertinent to z. Where prominence is irrelevant, NaN or zero values.
Each summit is condensed to a single value, although this may fail in some cases.
"""
function find_prominence(z, fofo; dic_neighbour = Dict{Symbol, String}())
    @assert isdir(fofo)
    # Retrieve boundary conditions.
    # This is a tuple of four vectors (from file, or zero-filled with the correct length)
    bcond = read_boundary_condition(fofo, size(z, 1), size(z, 2))
    @debug "    Find max tree "
    # We reduce the elevation resolution to whole meters in order to avoid oversegmentation,
    # which might lead to many uninteresting 'twin peaks' and heavier calculations. Many peaks are actually
    # split anyway.
    maxtree = MaxTree(round.(z))
    # Find the most important indices.
    @debug "    Find summits "
    summit_indices = distinct_summit_indices(z, maxtree)
    # Calculate mea
    @debug "    Find maximum elevation above"
    mea = maximum_elevation_above(z, bcond; maxtree, summit_indices)
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
    R_internal = CartesianIndices((3 : size(z, 1) - 2, 3 : size(z, 2) - 2))
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
    write_prominence_to_csv(prom, z, f_I_to_utm, indices, σ_indices, ffnam)
    ---> Bool

Writes elevation, prominence, position and index.

Return 'true' if the file was saved.
Return 'false' if the file would have the same elevation and prominence as previously.
"""
function write_prominence_to_csv(vpr, vz, vutm, indices, σ, ffnam)
    # Sheet indices
    vi = map(I -> I.I, indices)
    # We want to sort by
    #  1) prominence 'obscure' or  'prominent'
    #  2) elevation
    promlev_prom = get_config_value("Markers", "Prominence level [m], prominent summit", Int; nothing_if_not_found = false)
    sortval = [z + (p > promlev_prom ? 10000 : 0) for (z, p) in zip(vz, vpr)]
    # Order by elevation
    order = sortperm(sortval; rev = true)
    #
    # Prior to saving, look for changes compared to existing file.
    #
    if any_prominence_change_from_existing_csv_file(ffnam, vz[order], vpr[order])
        @debug "    Saving $ffnam"
        # Nest and prepare for output.
        # NOTE: The column order ought to remain the same as in
        # `add_names_to_csv`.
        vectors = [vz[order], vpr[order], vutm[order], vi[order], σ[order]]
        headers = ["Elevation_m", "Prominence_m", "Utm", "Sheet_index", "Stress"]
        widths = [20, 20, 20, 20, 20]
        write_vectors_to_csv(ffnam, headers, vectors, widths)
        return true
    else
        @info "No changes to $ffnam, skipping save."
        return false
    end
end

function any_prominence_change_from_existing_csv_file(ffnam, vz, vpr)
    if isfile(ffnam)
        old_data = readdlm(ffnam, '\t')[2:end,:]
        if isempty(old_data)
            if isempty(vz)
                return false
            else
                return true
            end
        end
        if size(old_data, 1) == length(vz)
            old_z, old_pr = old_data[:, 1], old_data[:, 2]
            if Float32.(old_z) == vz
                if Float32.(old_pr) == vpr
                    return false
                end
            end
        end
    end
    true
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
    mark_at!(img, prominent_indices, si2, sy2)
    # Single-pixel-width symbols (like a hollow triangle)
    # does not show up very clearly. We double it to increase 'line thicness'
    mark_at!(img, obscure_indices, si1, sy1)
    mark_at!(img, obscure_indices .+ CartesianIndex((1, 1)), si1, sy1)
    mark_at!(img, obscure_indices .+ CartesianIndex((-1, 1)), si1, sy1)
end

"""
    maximum_elevation_above(z::Array;
        maxtree = MaxTree(z),
        summit_indices = distinct_summit_indices(z, maxtree))

    maximum_elevation_above(z::Array, bcond;
        maxtree = MaxTree(z),
        summit_indices = distinct_summit_indices(z, maxtree))
    ---> similar(z)

Fill every summit's elevation region with its peak value.
If boundary conditions are supplied, regions may be dominated by boundary conditions.

We recommed rounding off and supply the keyword `maxtree = MaxTree(round.(z))`.

Also see 'prominence'.
"""
function maximum_elevation_above(z; maxtree = MaxTree(z), summit_indices = distinct_summit_indices(z, maxtree))
    bcond = boundary_cond_zero(z)
    maximum_elevation_above(z, bcond; maxtree, summit_indices)
end
function maximum_elevation_above(z::T, bcond;
    maxtree = MaxTree(z),
    summit_indices = distinct_summit_indices(z, maxtree)) where T <: AbstractArray
    #
    # Check arguments
    eltype(z) <: Float32 || throw(ArgumentError("eltype(z) = $(eltype(z)) is not <: Float32"))
    typeof(bcond) <: NTuple{8, Vector{Float32}} || throw(ArgumentError("typeof(bcond) = $(typeof(bcond)) is not <: NTuple{8, Vector{<: Float32}}"))
    # This algorithm should be more effective if we start with the tallest summits.
    # Consider that summit_indices is a Set. Make a sorted vector.
    vsummit_indices = let v = collect(summit_indices)
        v[sortperm(z[v]; rev = true)]
    end
    @assert typeof(bcond) <: NTuple{8, Vector{Float32}} "...typeof(bcond) = $(typeof(bcond))"
    # Allocate result matrix
    mea = fill(eltype(z)(NaN), size(z))
    # This takes 'mea' from neighbouring sheets into account.
    # Let's capture these variables in a closure:
    f! = func_mea_contact!(maxtree, z, bcond)
    # Visit direct and indirect parents of summits with
    # the information about the tallness of its tallest
    # descendants. The recursion returns when
    #    - a root is reached
    #    - or we reach a node that already knows about a taller parent.
    # At return, we pick the next node in this loop.
    for i in vsummit_indices
        max_elevation = z[i]
        parent_i = maxtree.parentindices[i]
        # Call the recursive, in_place function.
        f!(mea, max_elevation, parent_i)
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

"""
    add_names_to_csv(ffnam)

Lookup online, keep names in a sheet unique. Save as last column.
"""
function add_names_to_csv(ffnam)
    isfile(ffnam) || throw(ArgumentError("Does not exist: $ffnam"))
    old_data = readdlm(ffnam, '\t')[2:end,:]
    if isempty(old_data)
        @debug "    No summits in file => no names to add."
        return String[]
    end
    if size(old_data, 2) > 5
        @debug "    Names already exist => keep those unchanged."
        return String.(old_data[:, 6])
    end
    # utm positions as a vector of strings like "3,233"
    vsutm = map(s -> replace(s, ' ' => "", '(' => "", ')' => ""), old_data[:, 3])
    #
    endpoint = "/punkt"
    names = String[]
    for sutm in vsutm
        params = Dict(:koordsys => 25833,
            :utkoordsys => 25833,
            :nord => split(sutm, ',')[2],
            :ost => split(sutm, ',')[1],
            :radius => 150,
            :filtrer => "navn.stedsnavn,navn.meterFraPunkt")
        jsono = get_stadnamn_data(endpoint, params)
        if isempty(jsono)
            # Something went wrong. Don't write anything to file,
            # but return whatever we got so far.
            return names
        end
        lv1 = jsono.navn
        if isempty(lv1)
            push!(names, "")
        else
            if length(lv1) == 1
                lv2 = get(lv1[1], :stedsnavn)
            else
                # Take the closest place name. This
                # may not be the best choice always,
                # because e.g. the region name or anything
                # on a higher level might be assigned a closer coordinate
                # for arbitrary reasons.
                # NOTE, we could probably improve on this by checking navneobjekttype == "Fjell"
                closest_index = argmin(get.(lv1, "meterFraPunkt"))
                lv2 = get(lv1[closest_index], :stedsnavn)
            end
            if isempty(lv2)
                throw(ErrorException("Unexpected name, empty lv2, utm $sutm"))
            else
                if length(lv2) == 1
                    name = lv2[1].skrivemåte
                    if name ∈ names
                        @debug "    Encountered duplicate name $name without alternatives, utm $sutm. Name set to \"\""
                        push!(names, "")
                    else
                        push!(names, name)
                    end
                else
                    candidates = get.(lv2, :skrivemåte)
                    unique_candidates = filter(n-> n ∉ names, candidates)
                    if length(unique_candidates) == 1
                        push!(names, first(unique_candidates))
                    else
                        name = first(unique_candidates)
                        @debug "    Picked name \"$name\", at $sutm as the first in list: $(unique_candidates)."
                        push!(names, first(unique_candidates))
                    end
                end
            end
        end
    end
    length(names) == length(vsutm) || throw(ErrorException("Some names were not assigned. length $(length(names)) $(length(vsutm))"))
    @debug "    Saving names: $ffnam"
    # NOTE: The column order ought to remain the same as in
    # `write_prominence_to_csv`.
    vz, vpr, vutm, vi, vstress = old_data[:, 1], old_data[:, 2], old_data[:, 3], old_data[:, 4], old_data[:, 5]
    vectors = [vz, vpr, vutm, vi, vstress, names]
    headers = ["Elevation_m", "Prominence_m", "Utm", "Sheet_index", "Stress", "Name"]
    widths = [20, 20, 20, 20, 20, 20]
    write_vectors_to_csv(ffnam, headers, vectors, widths)
    names
end
