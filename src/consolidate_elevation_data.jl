
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


function consolidate_data_in_folder_to_geoarray(fofo)
    if isfile(joinpath(fofo, CONSOLIDATED_FNAM))
        @debug "$out_file_name in $fofo already exists. Exiting `consolidate_data_in_folder_to_geoarray`."
        return true
    end
    # The folder fofo's name contains the geometry request made at høydedata.no or similar!
    # r, c is row and column for this SheetBuilder of the entire SheetMatrixBuilder
    # The last four digits are the UTM corners used in the request.
    r, c, min_x, min_y, max_x, max_y = parse_folder_name(fofo)
    w = max_x - min_x
    h = max_y - min_y
    # Files to consolidate
    fnas_source = tif_full_filenames_buried_in_folder(fofo)
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
    # Displaying this for feedback might allocate as much memory as what we're trying to do. Still, for debugging:
    # display(map(g_dest.A[:,:, 1]) do pix
    #    nor = max(0.0, pix) / 1500
    #    PNGFiles.Gray(nor)
    #end)
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

Sample larger values occupying the same geographical position from g_source to g_dest.

The first method first checks if there is any overlap. If there is not, this method is faster.

The crs field is ignored, unlike `GeoArrays.sample_values!`.
"""
function sample_values_larger_than_limit!(g_dest::GeoArray, fna_source::String)
    g_source = GeoArrays.read(fna_source)
    if is_source_relevant(g_dest, g_source, fna_source)
        sample_values_larger_than_limit!(g_dest, g_source)
    end
    g_dest
end

function sample_values_larger_than_limit!(g_dest::GeoArray, g_source::GeoArray)
    wo, ho, zo = size(g_dest)
    w, h = size(g_source)[1:2]
    # Function that translates from logical coordinates
    # in `g_dest` to logical coordinates in `g_source`.
    f = inv(g_source.f) ∘ g_dest.f
    for io in 1:wo, jo in 1:ho
        # WAS i, j = round.(Int, f((io, jo)))
        i, j = Int.(f((io, jo)))
        # Is this logical coordinate inside the edges of g_source?
        if (1 <= i <= w) && (1 <= j <= h)
            # Loop over bands
            for z in 1:zo
                # Some source files may have value zero where they are really missing.
                # Also, close to zero elevation, there is much noise from waves and 
                # also from elevation recalibration. This heuristic gets rid of most.
                if g_source[i, j, z] > 0.7
                    g_dest[io, jo, z] = g_source[i, j, z]
                end
            end
        end
    end
    g_dest
end

function is_source_relevant(g_dest::T, g_source::T, fna_source::String) where T <: GeoArray
    is_source_relevant(bbox(g_dest), bbox(g_source), fna_source)
end
function is_source_relevant(d::T, s::T, fna_source) where T <: @NamedTuple{min_x::Float64, min_y::Float64, max_x::Float64, max_y::Float64}
    if !bbox_overlap(d, s)
        @debug "\t$(splitdir(fna_source)[end]) does not overlap $d"
        return false
    else
        return true
    end
end

