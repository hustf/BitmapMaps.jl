function consolidate_elevation_data(shp)
    consolidated = consolidate_data_in_folder_to_geoarray(full_folder_path(shp))
    @assert consolidated isa Bool
    if ! consolidated
        @info "Could not make consolidated .tif for sheet  with folder path $(shp.pthsh). Download and unzip .tif files? Exiting."
    end
    return consolidated
end


function consolidate_data_in_folder_to_geoarray(fofo)
    if isfile(joinpath(fofo, CONSOLIDATED_FNAM))
        @debug "$out_file_name in $fofo already exists. Exiting `consolidate_data_in_folder_to_geoarray`."
        return true
    end
    # The folder fofo's name contains the geometry request made at høydedata.no or similar!
    # r, c is row and column for this SheetBuilder of the entire SheetMatrixBuilder
    # The last four digits are the UTM corners used in the request.
    r, c, min_x, min_y, max_x, max_y = parse_folder_name(fofo)
    cbox = (;min_x, min_y, max_x, max_y)
    w = max_x - min_x
    h = max_y - min_y
    #
    fnas_source = tif_full_filenames_buried_in_folder(fofo)
    if length(fnas_source) == 0
        @info "No .tif files in $fofo to consolidate. Exiting."
        return false
    end
    g_dest = let
       # crs = GeoArrays.read(first(fnas_source)).crs
        A = zeros(Float32, w, h, 1)
        f = GeoArrays.AffineMap([1.0 0.0; 0.0 -1.0], 1.0 .* [min_x, max_y])
        GeoArray(A, f) # Note that when we droped defining crs here, we got  ERROR: AssertionError: 43200.5 down the line.
#        GeoArray(zeros(Float32, w, h - 1),  # TODO make very sure we want h - 1 and not h - 1
#            GeoArrays.AffineMap([1.0 0.0; 0.0 -1.0], 1.0 .* [min_x, max_y]), crs) # Note that when we droped defining crs here, we got  ERROR: AssertionError: 43200.5 down the line.
    end
    GeoArrays.bbox!(g_dest, (;min_x = min_x -1, min_y = min_y - 1, max_x, max_y))
    ii_min_x, ii_min_y = indices(g_dest, (cbox.min_x, cbox.min_y)).I
    ii_max_x, ii_max_y = indices(g_dest, (cbox.max_x, cbox.max_y)).I
    @show ii_min_x ii_min_y
    @show ii_max_x ii_max_y
    copy_sources_into_destination!(g_dest, min_x, max_y, max_x, min_y, fnas_source)
    # Displaying this for feedback might allocate as much memory as what we're trying to do. Still, for debugging:
    display(map(g_dest.A[:,:, 1]) do pix
        nor = max(0.0, pix) / 1500
        PNGFiles.Gray(nor)
    end)
    # Write to consolidated file
    GeoArrays.write(joinpath(fofo, CONSOLIDATED_FNAM), g_dest)
    return true
end

function copy_sources_into_destination!(g_dest, min_x, max_y, max_x, min_y, fnas_source)
    # Copy data file by file into g_dest
    for fna_source in fnas_source
        copy_source_into_destination!(g_dest, min_x, max_y, max_x, min_y, fna_source)
    end
    g_dest # by convention
end


function copy_source_into_destination!(g_dest, min_x, max_y, max_x, min_y, fna_source)
    #@show fna_source min_x max_x min_y max_y 
    #throw("sto")
    printstyled("\tPatching in $(splitdir(fna_source)[end])\n", color = :light_black)
    cbox = (;min_x, min_y, max_x, max_y)
    g_source = GeoArrays.read(fna_source)
    if ! is_geoarray_relevant(g_source,   cbox, fna_source)
        @debug "$fna_source is not relevant to $(cbox)"
    end
    # Find source corner index and geographical coordinate
    i1_source, j1_source = 1, 1
    i2_source, j2_source, _ = size(g_source.A)
    nw_source = coords(g_source, (i1_source, j1_source))
    # We check a lot because we don't know GeoArray very well...
    @assert nw_source[1] >= min_x "$(nw_source[1])"
    @assert nw_source[2] <= max_y + 0.5
    se_source = coords(g_source, (i2_source, j2_source))
    @assert se_source[1] <= max_x + 0.5
    @assert se_source[2] >= min_y
    @assert nw_source[1] <= se_source[1]
    @assert nw_source[2] >= se_source[2]
    # Find destination corner indices for these coordinates
    i1_dest, j1_dest = Tuple(indices(g_dest, nw_source)) # TODO add Vertex
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
    printstyled("\t\tCopying $(Int(round(length(v) / 1e6))) million points from $(splitdir(fna_source)[end]) \n\t at $(Dates.now())\n", color=:light_black)
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
    g_dest # by convention
end

"'boxwidth' does not include the width of the last little box at the end"
boxwidth(geoa) = boxwidth(bbox(geoa))
boxwidth(cbox::NamedTuple) = cbox.max_x - cbox.min_x
"'boxheight' does not include the height of the last little box at the end.
Assumes that northing is up."
boxheight(cbox::NamedTuple) = cbox.max_y - cbox.min_y
boxheight(geoa) = boxheight(bbox(geoa))



function is_geoarray_relevant(ga,   cbox, fna_source)
    if !bbox_overlap(bbox(ga), cbox)
        @debug "\t", splitdir(fna_source)[end], " does not overlap ", cbox
        return false
    elseif maximum(ga) <= 0.99
        @debug "\t", splitdir(fna_source)[end], " does not have elevations > 0.99"
        return false
    else
        #check_geoarray_details(ga; finam = splitdir(fna_source)[end])
        return true
    end
end

#= Dead 

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

=#