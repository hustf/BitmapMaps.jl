# Step in pipeline.
#
# Unlike most other pipeline functions, this takes captured arguments in addition to a SheetBuilder.
#
# We can't find summit prominence correctly without involving all sheets, since the taller summit may lie there.
# This step condenses the data from one sheet to a graph format, with which we can later combine sheets data.
#
# We output temporary Markers.png, Summits.csv and *.z files to this folder, and update the regional graph *.z file
# in the folder above.
#

"""
    summits_on_sheet(sb::SheetBuilder, ffna_graph, f_sides_with_border)
    summits_on_sheet(fofo, cell_iter, cell2utm, f_I_to_utm,
        prom_levels, symbols, symbol_sizes, min_stress,
        border_sides, ffna_graph)
    --> Bool

- Mark symbols in a temporary image MARKERS_FNAM.
- Summarize in a temporary .csv file for potential name tagging: SUMMITS_FNAM.
- Find the maxtree graph for every pixel, then extract data relevant to summits
- Merge the sheet's data with a regional graph, stored with file name `ffna_graph`, which finally
  will hold all the data required for finding summit prominence.
- Simplify the regional graph whenever its geographical outline is rectangular
"""
function summits_on_sheet(sb::SheetBuilder, ffna_graph, f_sides_with_border)
    # Harvest defaults
    min_prominence = get_config_value("Markers", "Prominence level [m], obscure summit", Int)
    min_stress = get_config_value("Markers", "Minimum stress level", Int)
    # 'Border sides' to this sheet actually depends on the sheet matrix builder. We have captured it in the closure:
    border_sides = f_sides_with_border(sb)
    # Arguments extracted from builders,, now do the work
    summits_on_sheet(full_folder_path(sb), sb.cell_iter, cell_to_utm_factor(sb), sb.f_I_to_utm, min_prominence, min_stress, border_sides, ffna_graph)
end
function summits_on_sheet(fofo, cell_iter, cell2utm, f_I_to_utm, min_prominence, min_stress, border_sides, ffna_graph)
    # Early exit:
    if ! isfile(joinpath(fofo, CONSOLIDATED_FNAM))
        @debug "    $CONSOLIDATED_FNAM in $fofo\n           does not exist. Exiting `summits_on_sheet`"
        return false
    end
    ffna_sum = joinpath(fofo, SUMMITS_FNAM)
    ffna_loc_graph = joinpath(fofo, splitdir(ffna_graph)[2])
    # Reuse earlier work if files exist.
    if ! isfile(ffna_sum)
        if isfile(ffna_loc_graph)
            @debug "    Creating first-pass $SUMMITS_FNAM  and replacing local elevation graph"
        else
            @debug "    Creating first-pass $SUMMITS_FNAM  and local elevation graph"
        end
        vI, maxtree, z = condense_summit_candidates_data(fofo, cell_iter, cell2utm, f_I_to_utm, min_stress, border_sides)
    elseif ! isfile(ffna_loc_graph)
        @debug "    Re-using $SUMMITS_FNAM  to build local elevation graph"
        # We shall add to regional elevation graph.
        # We can re-use the sheet indices from summits file:
        vI = read_indices_from_column(ffna_sum, 4)
        # Also add sheet borders to vI, regardless of if they are not included in the summits file.
        append_indices_of_borders!(vI, CartesianIndices(cell_iter), border_sides)
        # Load z for every cell (not the full source).
        z = elevation_at_output(fofo, cell_iter, cell2utm)
        # Limit z to minimum (will be rounded to 0)
        map!(ζ -> max(0.5f0, ζ), z, copy(z))
        # Calculating MaxTree once more is reasonably fast.
        maxtree = MaxTree(round.(z))
    else
        # Check if the regional graph is updated.
        if isfile(ffna_graph)
            gr = get_graph(ffna_graph)
            # Read the graph metadata with this syntax: gr[]
            already_merged = gr[]
            if ! isempty(already_merged)
                # Find out what needs to be done
                internal_borders = bbox_internal(cell_iter, f_I_to_utm)
                if internal_borders ∈ already_merged
                    # Late exit, but make sure Markers.png is there
                    # If Markers.png does not exist, make it, using squares to indicate temporary state.
                    if ! isfile(joinpath(fofo, MARKERS_FNAM))
                        @debug "    $MARKERS_FNAM is missing. Updating preliminary version from $SUMMITS_FNAM"
                        vI = read_indices_from_column(ffna_sum, 4)
                        update_marker_file(fofo, vI, cell_iter)
                    end
                    # Late exit, knowing that
                    # files in the folder are as nominal after this step.
                    @debug "    $SUMMITS_FNAM and local *.z file exists and was already merged with regional elevation graph. Exiting"
                    return true
                else
                    @debug "    $SUMMITS_FNAM and local *.z will be merged with existing regional elevation graph"
                end
            else
                throw(ErrorException("$(ffna_graph) exists, but contains no graph metadata. Consider deleting it, then start over."))
            end
        else
            @debug "    $SUMMITS_FNAM and local *.z files already exists"
        end
        gl = get_graph(ffna_loc_graph)
    end
    if ! @isdefined gl
         @debug "    Building local graph from maxtree"
         # Add summits and relevant border to the regional elevation graph.
        # Many of the summit prominences can be determined locally,
        # and those with low prominece are filtered out in this first pass.
        # When more regions are added, we can filter out more low-prominence summits.
        gl = build_and_save_sheet_graph(ffna_loc_graph, maxtree, z, vI, f_I_to_utm, border_sides, min_prominence)
    end
    # Use the pruned local graph to reduce the number of summits in Summits.csv
    # Store the reduced vector of cell indices, vIr.
    vIr = update_summits_file(ffna_sum, gl, f_I_to_utm, cell_iter)
    # Update Markers.png, using squares to indicate the temporary state.
    # A large number of stones and trees will be close to the border in the graph
    update_marker_file(fofo, vIr, cell_iter)
    #
    # Merge the sheet local graph into the region graph
    #
    # Get the (larger) regional graph, if it was not already loaded to check metadata.
    if ! @isdefined gr
        gr = get_graph(ffna_graph)
    end
    # Merge
    @debug "    Merging local into regional graph, starting with $(nv(gr)) regional vertices"
    merge_into!(gr, gl, cell_iter, f_I_to_utm, border_sides)
    # If `border_sides` does not include :n, this sheet is the topmost in its column.
    # Which also means that the regional graph  is rectangular.
    # And if the regional graph is rectangular, we can apply `reduce_and_prune!`,
    # without harming the external connections. `reduce_and_prune!` will remove 'diamonds'
    # on the sheet borders which are not not external to the region.
    if :n ∉ border_sides || :w ∉ border_sides
        # When more sheets are added to the regional graph, we can determine
        # the prominence for more of the summits. Hence, we can cull some more
        # even if the graph is not complete. But we only do that when the regional
        # graph is in a state where it is rectangular.
        @debug "    Regional outline is rectangular.  `reduce_and_prune!`"
        # We won't prune vertices on the sides that are yet to be connected.
        # Those sides will be protected at this point.
        if :w ∉ border_sides && :n ∈ border_sides
            protected_sides = intersect(border_sides, [:e, :n])
        else
            protected_sides = intersect(border_sides, [:e])
        end
        limiting_bbox = grow_bbox_on_sides_excepting(bbox_internal(gr), protected_sides)
        reduce_and_prune!(gr, limiting_bbox, min_prominence)
    end
    # Save the larger regional graph (in MGFormat by default)
    savegraph(ffna_graph, gr)
    true
end

"""
    update_summits_file(ffna_sum, g, f_I_to_utm, cell_iter)
    --> Vector{CartesianIndex{2}}

Eliminate lines in the summits file `ffna_sum` by keeping those lines that are
leaves in the graph `g` (but do not lie on the border of the sheet).

Return cell indices for the kept summits.
"""
function update_summits_file(ffna_sum, g, f_I_to_utm, cell_iter)
    # Find coordinates of summits on this sheet, except "summits" on borders
    # with other sheets.
    crop_west, crop_north = f_I_to_utm(CartesianIndex(2, 2))
    crop_east, crop_south,  = f_I_to_utm(CartesianIndex(cell_iter[end].I .- (1, 1)))
    set_utm = Set(filter!(leaf_utms(g)) do (easting, northing)
        easting >= crop_west &&
        easting <= crop_east &&
        northing >= crop_south &&
        northing <= crop_north
    end)
    # Read the .csv file to be filtered.
    headers = string.(split(readline(ffna_sum), '\t'))
    summits_data = readdlm(ffna_sum, '\t')[2:end,:]
    @assert size(summits_data, 2) >= 5
    # Early exit.
    if length(headers) == 6
        @debug "    Skip updating summits file from local graph since names are already added"
        # Return the cell indices of the summits already in file.
        return read_indices_from_column(ffna_sum, 4)
    end
    pruned_summits_data = filter(eachrow(summits_data)) do (z, pr, utm_string, I, σ)
        utm = parse_as_int_tuple(utm_string)
        utm ∈ set_utm
    end
    # For storage in a modified .csv file, we'll nest by columns instead of by rows.
    # Initialize the new structure, then populate it.
    if ! isempty(pruned_summits_data)
        vectors = [Vector{typeof(pruned_summits_data[1][i])}(undef, length(pruned_summits_data)) for i in 1:5]
        # Populate the vectors.
        for i in 1:length(pruned_summits_data)
            for j in 1:5
                vectors[j][i] = pruned_summits_data[i][j]
            end
        end
    else
        # Oh, there are no summits on this sheet. Probably ocean?
        vectors = [Float64[], Float64[], SubString{String}[], SubString{String}[], Float64[]]
    end
    widths = repeat([20], length(vectors))
    write_vectors_to_csv(ffna_sum, headers, vectors, widths)
    # Return the cell indices of the kept summits.
    CartesianIndex.(parse_as_int_tuple.(vectors[4]))
end


"""
    update_marker_file(fofo, vI, cell_iter)
    update_marker_file(fofo, vI, cell_iter, prom_levels, symbols, symbol_sizes)

Use the last method when summit prominence is actually known exactly, which we can't
before all sheets have been processed once.
"""
function update_marker_file(fofo, vI, cell_iter)
    # Save a temporary markers image. This will be changed later, after regional prominence calculation.
    # We don't know the last three arguments yet.
    update_marker_file(fofo, vI, cell_iter, [0, 1], ["in_square", "in_square"], [11, 11])
end
function update_marker_file(fofo, vI, cell_iter, prom_levels, symbols, symbol_sizes)
    mffna = joinpath(fofo, MARKERS_FNAM)
    # We use squares to make the temporary nature of this image evident.
    # We don't know any prominence values yet, so all are simply zero.
    bwres = draw_summit_marks(map(x -> 0.0f0, cell_iter), vI, prom_levels, symbols, symbol_sizes)
    # Feedback
    display_if_vscode(bwres)
    @debug "    Saving $mffna"
    # Make a transparent color image, save it.
    save_png_with_phys(mffna, map(bwres) do pix
        pix == true && return RGBA{N0f8}(0., 0, 0, 1)
        RGBA{N0f8}(0., 0, 0, 0)
    end)
end

"""
    condense_summit_candidates_data(fofo, cell_iter, cell2utm, f_I_to_utm, min_stress, border_sides)
    --> vI, maxtree, z

Finding summit prominence is a two-stage process, since calculating actual prominence requires info from, in principle, all sheets.
This is the first step:

- Identify candidates for summits from elevation (will be filtered down in later steps)
- Summarize in .csv file: SUMMITS_FNAM.
#- Mark symbols in output image:  MARKERS_FNAM. TODO: Take this from graph, not from Summits.csv
- Update relevant values in regional (sheet matrix) graph: ffna_graph

The relevant indices for the graph is: summit candidates and sheet borders where they are adjacent to other sheets.
"""
function condense_summit_candidates_data(fofo, cell_iter, cell2utm, f_I_to_utm, min_stress, border_sides)
    ny, nx = size(cell_iter)
    # We need the full source matrix in this function. It's size is size(cell_iter) * cell2utm^2,
    # so performance is low with this function.
    zf = elevation_full_res(fofo)
    # Limit z to minimum (will be rounded to 0)
    map!(ζ -> max(0.5f0, ζ), zf, copy(zf))
    # Source index into the transposed source
    si = CartesianIndices((1:cell2utm:(ny  * cell2utm), 1:cell2utm:(nx * cell2utm)))
    #
    # We keep the full value precision here, so that we
    # are able to pick the very topmost cell.
    z = zf[si]
    R = CartesianIndices(z)
    # Make the maxtree, which is an effective elevation graph format for images.
    # We reduce the elevation precision to whole meters here. This is in order to avoid oversegmentation,
    # which might lead to many uninteresting 'twin peaks' and heavier calculations. Many peaks are actually
    # split anyway.
    maxtree = MaxTree(round.(z))
    # Find a rather large set of summit indices(roughly 1 in 100 out of all z indices).
    # Each summit is a linear index into z.
    # This set does not contain summits in immediate proximity to each other.
    summit_i_set = distinct_summit_indices(z, maxtree)
    # Reduce the set a little:
    #   1) Drop summit elevation < 200 m, hardcoded. This removes a lot of power-line "summits", waves and other artifacts.
    #      (we now may have 1 in 150 out of all R)
    filter!(i -> z[i] > 200, summit_i_set)
    #
    #
    #   2) Remove some more false peaks based on user-defined min_stress,
    #      i.e. very large negative amplitude Hessian components.
    #      This removes artifacts, most power lines, and large conifers.
    #      Most large conifer trees were already filtered out by the elevation filter above.
    #
    #  Tensile stress is negative on top of summits, positive in valleys.
    #  We use the full resolution for calculating stress, so results won't be affected by cell2utm.
    #
    #  Unsmoothed mean stress, downsampled from the full resolution source data
    σ = σm(zf)[si]
    filter!(i -> 0 > σ[i] > min_stress, summit_i_set)
    # Make an ordered vector of CartesianIndex from the large summit candidates list.
    vI = R[collect(summit_i_set)]
    # We have not as of yet an idea of the peak prominence, which may well be a metre only.
    # For calculation of prominence (often influenced by neigbouring sheets), we need
    # vertices on the borders which has a neighbour. These will be connected to neighbouring
    # sheet while the local graph is merged into the regional graph.
    append_indices_of_borders!(vI, R, border_sides)
    #
    # The stress result is stored in a corresponding vector. This will be useful for results inspection
    # after further filtering.
    vσ = [σ[I] for I in vI]
    #
    # Store the results. Most of these summit candidates are probably trees, and will be discarded later
    # due to low prominence. However, we want to avoid calculating the mean stress again since it takes ~15 seconds.
    waschanged = write_prominence_to_csv(zeros(length(vI)), round.(z[vI]), map(f_I_to_utm, vI), vI, vσ, joinpath(fofo, SUMMITS_FNAM))
    vI, maxtree, z
end



function draw_summit_marks(prominence, indices, prom_levels, symbols, symbol_sizes)
    length(prom_levels) == length(symbols) == length(symbol_sizes) == 2 || throw(ArgumentError("Length not 2"))
    #
    @debug "    Draw summit marks $symbols"
    # Equilateral triangles are drawn with a horizontal line at bottom.
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
    isempty(prominent_indices) || mark_at!(img, prominent_indices, si2, sy2)
    if ! isempty(obscure_indices)
        # Single-pixel-width symbols (like a hollow triangle)
        # do not show up very clearly. We double it to increase 'line thickness'
        mark_at!(img, obscure_indices, si1, sy1)
        mark_at!(img, obscure_indices .+ CartesianIndex((1, 1)), si1, sy1)
        mark_at!(img, obscure_indices .+ CartesianIndex((-1, 1)), si1, sy1)
    end
    img
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

    Case A returns Set(5), cartesian index [2, 2]:

```
julia> distinct_summit_indices(A, MaxTree(A))
Set{Int64} with 1 element:
  5
```

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
    # For each family, add the index of its tallest member,
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
    promlev_prom = get_config_value("Markers", "Prominence level [m], prominent summit", Int)
    sortval = [z + (p > promlev_prom ? 10000 : 0) for (z, p) in zip(vz, vpr)]
    # Order
    order = sortperm(sortval; rev = true)
    #
    # Prior to saving, look for changes compared to existing file.
    #
    if any_prominence_change_from_existing_csv_file(ffnam, vz[order], vpr[order])
        @debug "    Saving $ffnam"
        # Nest and prepare for output.
        # TODO: Use `write_named_tuple_to_csv`
        vectors = [vz[order], vpr[order], vutm[order], vi[order], σ[order]]
        headers = ["Elevation_m", "Prominence_m", "Utm", "Cell_index", "Stress"]
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



"""
    append_indices_of_borders!(vI, R, border_sides)

Append indices of borders, return unique indices.

vI is a vector of CartesianIndices in R.

border_sides is e.g. [:n, :s, :e], indicating North (up), etc.

Corners are left out so as not to make the regional elevation graph overly connected.
"""
function append_indices_of_borders!(vI, R, border_sides)
    foreach(border_sides) do border_side
        fullside = indices_of_border(R, border_side)
        # Don't include any end
        internalside = fullside[2:end - 1]
        append!(vI, internalside)
    end
    vI
end

"""
    indices_of_border(R, border_side)
    --->  Vector{CartesianIndex{2}}

R is CartesianIndices, e.g.
```
julia> R = CartesianIndices((5500, 3820));

julia> Rutm = CartesianIndices((25785:2:33423, 6927028:-2:6916030))
```

Also see `pairs_of_border_neighbors`.
"""
function indices_of_border(R, border_side)
    if :n == border_side
        R[1, :]
    elseif :s == border_side
        R[end, :]
    elseif :w == border_side
        R[:, 1]
    elseif :e == border_side
        R[:, end]
    else
        throw(ArgumentError("border_side is $border_side, not :n, :s, :e:, :w"))
    end
end

"""
    read_indices_from_column(ffna_sum, column_no)
    ---> Vector{CartesianIndex{2}}
"""
function read_indices_from_column(ffna_sum, column_no)
    @assert column_no ∈ [3, 4]
    summits_index_strings = readdlm(ffna_sum, '\t')[2:end, column_no]
    map(summits_index_strings) do index_string
        CartesianIndex(parse_as_int_tuple(index_string))
    end
end