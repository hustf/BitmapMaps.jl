# Step in pipeline
"""
    consolidate_elevation_data(sb) ---> Bool
"""
function consolidate_elevation_data(sb)
    consolidated = consolidate_data_in_folder_to_geoarray(full_folder_path(sb))
    @assert consolidated isa Bool
    if ! consolidated
        @info "Could not make consolidated .tif for sheet  with folder path $(sb.pthsh). Download and unzip .tif files? Exiting."
    end
    return consolidated
end

"""
    consolidate_data_in_folder_to_geoarray(fofo)
    ---> Bool
"""
function consolidate_data_in_folder_to_geoarray(fofo)
    if isfile(joinpath(fofo, CONSOLIDATED_FNAM))
        @debug "    $CONSOLIDATED_FNAM in $fofo already exists. Exiting `consolidate_data_in_folder_to_geoarray`"
        return true
    end
    # The folder fofo's name contains the geometry request made at høydedata.no or similar!
    # r, c is row and column for this SheetBuilder of the entire SheetMatrixBuilder
    # The last four digits are the UTM corners used in the request.
    r, c, min_x, min_y, max_x, max_y = parse_folder_name(fofo)
    w = max_x - min_x
    h = max_y - min_y
    # Files to consolidate
    fnas_source = filter(tif_full_filenames_buried_in_folder(fofo)) do ffna
        splitpath(ffna)[end] !== CONSOLIDATED_FNAM
    end
    if length(fnas_source) == 0
        @info "No .tif files in $fofo to consolidate. Exiting."
        return false
    end
    g_dest = let
        A = zeros(Float32, w, h, 1)
        f = GeoArrays.AffineMap([1.0 0.0; 0.0 -1.0], 1.0 .* [min_x, max_y])
        GeoArray(A, f)
    end
    copy_sources_into_destination!(g_dest, fnas_source)
    if sum(g_dest.A) == 0
        @info "No relevant data to consolidate. Exiting."
        return false
    end
    # Write to consolidated file
    GeoArrays.write(joinpath(fofo, CONSOLIDATED_FNAM), g_dest)
    return true
end

function copy_sources_into_destination!(g_dest, fnas_source)
    # Copy data file by file into g_dest
    for fna_source in fnas_source
        sample_values_larger_than_limit!(g_dest, fna_source)
    end
    g_dest # by convention
end

"""
    sample_values_larger_than_limit!(g_dest::GeoArray, fna_source::String)
    sample_values_larger_than_limit!(g_dest::GeoArray, g_source::GeoArray)

Sample values occupying the same geographical position from g_source to g_dest.

The crs field is ignored, unlike `GeoArrays.sample_values!`.
"""
function sample_values_larger_than_limit!(g_dest::GeoArray, fna_source::String)
    sample_values_larger_than_limit!(g_dest, readclose(fna_source))
end

function sample_values_larger_than_limit!(g_dest::GeoArray, g_source::GeoArray)
    wo, ho, zo = size(g_dest)
    w, h = size(g_source)[1:2]
    # Function that translates from logical coordinates
    # in `g_dest` to logical coordinates in `g_source`.
    f = inv(g_source.f) ∘ g_dest.f
    for io in 1:wo, jo in 1:ho
        i, j = Int.(f((io, jo)))
        # Is this logical coordinate inside the edges of g_source?
        if (1 <= i <= w) && (1 <= j <= h)
            # Loop over bands
            for z in 1:zo
                # Some source files has value zero where data is really missing.
                # Also, close to zero elevation, there is much noise from waves and
                # also from elevation recalibration. This heuristic gets rid of most.
                val = g_source[i, j, z]
                if ! ismissing(val) && val > 0.7
                    g_dest[io, jo, z] = g_source[i, j, z]
                end
            end
        end
    end
    g_dest
end

function is_source_relevant(g_dest::GeoArray, fna_source::String)
    # Destination bounding box is determined including zero-valued cells.
    bb_dest = bbox(g_dest)
    min_x = Int(bb_dest.min_x)
    min_y = Int(bb_dest.min_y)
    max_x = Int(bb_dest.max_x)
    max_y = Int(bb_dest.max_y)
    bbd = (;min_x, min_y, max_x, max_y)
    # Do boxes overlap (adjacent is not overlap)?
    is_source_relevant(bbd, fna_source)
end

function is_source_relevant(bb_dest, fna_source::String)
    # Source bounding box neglects zero-padded (empty) cells.
    bb_source = nonzero_raster_rect(fna_source)
    # Do boxes overlap (adjacent is not overlap)?
    bbox_external_overlap(bb_dest, bb_source)
end