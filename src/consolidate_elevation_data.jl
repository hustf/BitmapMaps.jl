# Step in pipeline.
# Pulls data from relevant .tif geodata files into a single
# file stored in the sheet builder's directory.
# Reading from a consolidated file is faster because
# we can assume the numeric type (and no Missing) values,
# and because we won't load data external to the sheet.

"""
    consolidate_elevation_data(sb::SheetBuilder)
    ---> Bool

`sb` is a SheetBuilder, which is part of a collection, a SheetMatrixBuilder.

This calls `consolidate_local_data_to_geoarray_in_folder(sb.pthsh)`.
See that function regarding where to place input .tif files.
"""
function consolidate_elevation_data(sb)
    consolidated = consolidate_local_data_to_geoarray_in_folder(full_folder_path(sb); include_parent_folder = true)
    @assert consolidated isa Bool
    if ! consolidated
        @info "Could not make $CONSOLIDATED_FNAM for sheet  with folder path $(sb.pthsh). Download and unzip .tif files? Exiting."
    end
    return consolidated
end

"""
    consolidate_local_data_to_geoarray_in_folder(fofo; include_parent_folder = false)
    ---> Bool

The folder name 'fofo' is interpreted as the external bounding box.
Data which fit into the box is collected from any .tif files in the folder itself.
Provided that 'include_parent_folder = true', the "sheet matrix folder" is also searched.

'fofo' is intended as a working directory for one sheet in the total BitmapMap, specified
by a SheetBuilder, one of normally several in a SheetMatrixBuilder.

# Example of folder structure

If 'fofo1' and 'fofo2' is contained in 'fo', it is good practice to name 'fo' similarly.  e.g.:

´´´
fo:    homedir()/BitmapMaps/myproj 47675 6929520 50858 6938686
fofo1:                            \1 1  47675 6929520  50858 6934103
fofo2:                            \2 1  47675 6934103  50858 6938686
´´´

If the above structure is used, and the collected .tif data files cover a larger area than single sheets,
drop the data files in 'fo' and let this function pull data from 'fo' into each 'fofo/Consolidated.tif`.

The same procedure would regardless of the collected data file's extent (but the consolidation would take
more than if data files were put straight in the correct 'fofo'.).

Also see `copy_relevant_tifs_to_folder`.
"""
function consolidate_local_data_to_geoarray_in_folder(fofo; include_parent_folder = false)
    if isfile(joinpath(fofo, CONSOLIDATED_FNAM))
        @debug "    $CONSOLIDATED_FNAM in $fofo \n           already exists. Exiting `consolidate_local_data_to_geoarray_in_folder`"
        return true
    end
    # The folder fofo's name contains the geometry request made at høydedata.no or similar!
    # r, c is row and column for this SheetBuilder of the entire SheetMatrixBuilder
    # The last four digits are the UTM corners used in the request.
    r, c, min_x, min_y, max_x, max_y = parse_folder_name(fofo)
    w = max_x - min_x
    h = max_y - min_y
    # Files to consolidate
    alltifs = tif_full_filenames_buried_in_folder(fofo)
    if include_parent_folder
        # We may previously have examined the source files. This is time consuming, so let's update it from file.
        read_TIFDIC(abspath(joinpath(fofo, "..")))
        append!(alltifs, tif_full_filenames_in_parent_folder(fofo))
        save_TIFDIC(abspath(joinpath(fofo, "..")))
    end
    fnas_source = filter(alltifs) do ffna
        splitpath(ffna)[end] !== CONSOLIDATED_FNAM
    end
    if length(fnas_source) == 0
        @info "No .tif files in $fofo to consolidate. Exiting."
        return false
    end
    g_dest = let
        A = zeros(Float32, w, h)
        linear = SA[1.0 0.0; 0.0 -1.0]
        translation = SA[Float64(min_x), Float64(max_y)]
        f = AffineMap(linear, translation)
        GeoArray(A, f)
    end
    copy_sources_into_destination!(g_dest, filter(fna -> is_source_relevant(g_dest, fna), fnas_source))
    if sum(g_dest.A) == 0
        allow_empty_sheets = get_config_value("Behaviour when data is missing", "Fill with elevation zero (true or false)", String)
        if allow_empty_sheets == "false"
            @info "No relevant data to consolidate. \n      Consider changing section 'Behaviour when data is missing' in $(_get_fnam_but_dont_create_file()). Exiting."
            return false
        elseif allow_empty_sheets !== "true"
            throw(ArgumentError("Unexpected value of 'allow_empty_sheets': $(allow_empty_sheets)"))
        end
    end
    # Feedback
    display_if_vscode(permutedims(g_dest.A[:, :]))
    # Write to consolidated file
    GeoArrays.write(joinpath(fofo, CONSOLIDATED_FNAM), g_dest)
    return true
end

function copy_sources_into_destination!(g_dest, fnas_source)
    # Copy data file by file into g_dest
    for fna_source in fnas_source
        @debug "    Copy cells from $fna_source"
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
    g_source = readclose(fna_source)
    if cell_to_utm_factor(g_source) > 1
        @warn "cell_to_utm_factor($(fna_source)) > 1"
    end
    sample_values_larger_than_limit!(g_dest, g_source)
end

function sample_values_larger_than_limit!(g_dest::GeoArray, g_source::GeoArray)
    wo, ho = size(g_dest)
    w, h = size(g_source)
    # Function that translates from logical coordinates
    # in `g_dest` to logical coordinates in `g_source`.
    f = inv(g_source.f) ∘ g_dest.f
    for io in 1:wo, jo in 1:ho
        i, j = Int.(round.(f((io, jo))))
        # Is this logical coordinate inside the edges of g_source?
        if (1 <= i <= w) && (1 <= j <= h)
            # Some source files has value zero where data is really missing.
            # Also, close to zero elevation, there is much noise from waves and
            # also from elevation recalibration. This heuristic gets rid of most.
            val = g_source[i, j]
            if ! ismissing(val)
                g_dest[io, jo] = g_source[i, j]
            end
        end
    end
    g_dest
end

function is_source_relevant(g_dest::GeoArray, fna_source::String)
    # Destination bounding box is determined including zero-valued cells.
    bb_dest = bounds(bbox(g_dest))
    min_x, max_x = Int.(bb_dest.X)
    min_y, max_y = Int.(bb_dest.Y)
    bbd = (;min_x, min_y, max_x, max_y)
    # Do boxes overlap (adjacent is not overlap)?
    is_source_relevant(bbd, fna_source)
end

function is_source_relevant(bb_dest, fna_source::String)
    # Source bounding box neglects zero-padded (empty) cells.
    @debug "    Find nonzero bounding box from $fna_source"
    bb_source = nonzero_raster_rect(fna_source)
    # Do boxes overlap (adjacent is not overlap)?
    bbox_external_overlap(bb_dest, bb_source)
end