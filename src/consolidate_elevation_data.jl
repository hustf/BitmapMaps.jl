function consolidate_elevation_data(shp)
    consolidated = consolidate_data_in_folder_to_geoarray(full_folder_path(shp))
    @assert consolidated isa Bool
    if ! consolidated
        @info "Could not make consolidated .tif for sheet partitioning with folder path $(shp.pthsh). Download and unzip .tif files? Exiting."
    end
    return consolidated
end


function consolidate_data_in_folder_to_geoarray(fofo)
    if isfile(joinpath(fofo, CONSOLIDATED_FNAM))
        @debug "$out_file_name in $fofo already exists. Exiting `consolidate_data_in_folder_to_geoarray`."
        return true
    end
    # The folder fofo's name contains the geometry request made at høydedata.no or similar!
    # r, c is row and column for this SheetPartition of the entire BmPartition
    # The last four digits are the UTM corners used in the request.
    r, c, min_x, min_y, max_x, max_y = parse_folder_name(fofo)
    w = max_x - min_x
    h = max_y - min_y
    cbox = (;min_x, min_y, max_x, max_y)
    #
    fis = tif_full_filenames_buried_in_folder(fofo)
    if length(fis) == 0
        @info "No .tif files in $fofo to consolidate. Exiting."
        return false
    end
    # The read arrays may go outside requested area. We read and then crop each,
    # dropping those that concern depth under water
    relevant_geoas, relevant_tif_names = relevant_cropped_geoas_and_filenames(fis, cbox)
    g_dest = GeoArray(zeros(Float32, w, h - 1), 
        GeoArrays.AffineMap([1.0 0.0; 0.0 -1.0], 1.0 .* [min_x, max_y]),
        first(relevant_geoas).crs)
    check_geoarray_details(g_dest; min_x, max_y, max_x, min_y)
    copy_checked_data_files_to_destination!(g_dest, min_x, max_y, max_x, min_y, relevant_geoas, relevant_tif_names)
    # Displaying this for feedback might allocate as much memory as what we're trying to do. Still, for debugging:
    display(map(g_dest.A[:,:, 1]) do pix
        nor = max(0.0, pix) / 1500
        PNGFiles.Gray(nor)
    end)
    # Write to consolidated file
    GeoArrays.write(joinpath(fofo, CONSOLIDATED_FNAM), g_dest)
    return true
end

function copy_checked_data_files_to_destination!(g_dest, min_x, max_y, max_x, min_y, relevant_geoas, relevant_tif_names)
    # Now move the info from each relevant geoa source into g_dest.
    for (g_source, fi) in zip(relevant_geoas, relevant_tif_names)
        copy_checked_data_to_destination!(g_dest, min_x, max_y, max_x, min_y, g_source, fi)
    end
end


function copy_checked_data_to_destination!(g_dest, min_x, max_y, max_x, min_y, g_source, fi)
    @show fi min_x max_x min_y max_y 
   # throw("sto")
    printstyled("\tPatching in $(splitdir(fi)[end])\n", color = :light_black)
    # Find source corner index and geographical coordinate
    i1_source, j1_source = 1, 1
    i2_source, j2_source, _ = size(g_source.A)
    nw_source = coords(g_source, (i1_source, j1_source))
    # We check a lot because we don't know GeoArray very well...
    @assert nw_source[1] >= min_x
    @assert nw_source[2] <= max_y + 0.5
    se_source = coords(g_source, (i2_source, j2_source))
    @assert se_source[1] <= max_x + 0.5
    @assert se_source[2] >= min_y
    @assert nw_source[1] <= se_source[1]
    @assert nw_source[2] >= se_source[2]
    # Find destination corner indices for these coordinates
    i1_dest, j1_dest = Tuple(indices(g_dest, nw_source))
    i2_dest, j2_dest = Tuple(indices(g_dest, se_source))
    # Limit because we may miss a little. Don't know why, don't think it matters.
    i1_dest < 1 && @warn "i1_dest" i1_dest 
    j1_dest < 1 && @warn "j1_dest" j1_dest 
    i1_dest = max(i1_dest, 1)
    j1_dest = max(j1_dest, 1)
    i2_dest > size(g_dest, 1) && @warn "i2_dest" i2_dest 
    j2_dest < size(g_dest, 2) && @warn "j2_dest" j2_dest 
    i2_dest = min(i2_dest, size(g_dest, 1))
    j2_dest = min(j2_dest, size(g_dest, 2))
    # Maybe not the intended way to add this, but it works, and naive broadcasting (.+=) does not.
    v = @view g_dest[i1_dest:i2_dest, j1_dest:j2_dest, begin:end];
    printstyled("\t\tCopying $(Int(round(length(v) / 1e6))) million points from $(splitdir(fi)[end]) \n\t at $(Dates.now())\n", color=:light_black)
    for i in eachindex(v)
        # We do suspect there may be a row or column overlapping between source files. Which is fine 
        # with this high resolution. 
        if v[i] == 0.0
            # There is a good deal of areas close to zero.
            # That looks really bad on topographic reliefs.
            # We can hide some, accepts some, but discard these most common!
            if g_source[i] > 0.7
                v[i] = g_source[i]
            end
        end
    end
    printstyled("\t\tOk\n", color=:light_black)
    g_dest
end

"'boxwidth' does not include the width of the last little box at the end"
boxwidth(geoa) = (last(coords(geoa)) - first(coords(geoa)))[1]
boxwidth(cbox::NamedTuple) = cbox.max_x - cbox.min_x
"'boxheight' does not include the height of the last little box at the end.
Assumes that northing is up."
boxheight(cbox::NamedTuple) = cbox.max_y - cbox.min_y
boxheight(geoa) = boxheight(bbox(geoa))



function check_geoarray_details(g; min_x = 0, max_y = 0, max_x = 0, min_y = 0, finam = "")
    if ! isempty(finam)
        printstyled("\tChecking $finam    ", color = :light_black)
    end
    # Start corner indices -> Utm coordinates of north-west corner
    we, no = coords(g, (1, 1), Vertex())
    min_x !== 0 && @assert we == min_x
    max_y !== 0 && @assert no == max_y 
    # Utm coordinates of north-west corner -> Start corner indices 
    @assert Tuple(indices(g, (we, no), Vertex())) == (1,1)
    # Increasing easting -> increasing column index.
    # This also checks that the resolution is 1 utm <=> 1 index
    @assert Tuple(indices(g, (we + 1, no), Vertex())) == (2,1)
    # Decreasing northing -> increasing row index
    @assert Tuple(indices(g, (we, no - 1), Vertex())) == (1,2)
    # Any value in the matrix represents a range of coordinates - a square on the map.
    # When asking for the coordinate of these indices, we get the west-south corner
    # of that square.
    @assert g.f.translation == [we, no]
    w = boxwidth(g)
    h = boxheight(g)
    ie, je, _ = size(g)
    ea, so = coords(g, (ie, je), Vertex())
    @assert ea - we == w
    @assert no - so == h - 1 # Note how strange and wonderful.
    @assert Tuple(indices(g, (ea, so), Vertex())) == (ie, je)
end


function relevant_cropped_geoas_and_filenames(fis, cbox)
    relevant_geoas = GeoArray[]
    relevant_tif_names = String[]
    for fi in fis
        ga = GeoArrays.read(fi)
        if !bbox_overlap(bbox(ga), cbox)
            printstyled("\t", splitdir(fi)[end], " does not overlap ", cbox, "\n", color = :light_black)
            continue
        elseif maximum(ga) <= 0.99
            printstyled("\t", splitdir(fi)[end], " does not have elevations > 0.99 \n", color = :light_black)
            continue
        else
            check_geoarray_details(ga; finam = splitdir(fi)[end])
            println()
            gc = crop(ga, cbox)
            push!(relevant_geoas, gc)
            push!(relevant_tif_names, fi)
        end
    end
    relevant_geoas, relevant_tif_names
end
