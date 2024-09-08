# Step in pipeline.
# Identifying water surfaces is a time consuming operation.
# Results look good, but may benefit from manual touch-up of the output file.
# This might probably be speed up a good deal.
"""
    water_overlay(sb::SheetBuilder)
    water_overlay(fofo)
    ---> Bool

Creates water overlay from elevation data.
Output is an RGBA .png image file. Fully transparent outside water surfaces, for manual touch-up.
"""
function water_overlay(sb::SheetBuilder)
    lake_steepness_max = get_config_value("Water", "Lake steepness max", Float32)
    # Go ahead
    water_overlay(full_folder_path(sb), sb.cell_iter, cell_to_utm_factor(sb), lake_steepness_max)
end
function water_overlay(fofo, cell_iter, cell2utm, lake_steepness_max)
    if ! isfile(joinpath(fofo, CONSOLIDATED_FNAM))
        @debug "    $CONSOLIDATED_FNAM in $fofo\n           does not exist. Exiting `water_overlay`"
        return false
    end
    ffnam_csv_lakes = splitext(joinpath(fofo, WATER_FNAM))[1] * ".csv"
    if isfile(joinpath(fofo, WATER_FNAM)) && isfile(ffnam_csv_lakes) 
        @debug "    $WATER_FNAM in $fofo\n           already exists. Exiting `water_overlay`"
        return true
    end
    ffna = joinpath(fofo, CONSOLIDATED_FNAM)
    elevation = let
        g = readclose(ffna)
        eltype(g) == Float32 || throw(TypeError(:g, "unexpected .tif image eltype", GeoArrays.GeoArray{Float32, Array{Float32, 3}}, typeof(g)))
        # We're transposing the source data here, because
        # it makes it easier to reason about north, south, east west.
        transpose(g.A[:, :, 1])
    end
    # Find segments with lakes, and without. Whether a segment is a lake
    # determines segment_mean(..)
    segments = lake_segment(elevation, cell_iter, cell2utm, lake_steepness_max)
    # Save lakes posittions and elevations in a .csv file
    write_lake_to_csv(ffnam_csv_lakes , segments, elevation, cell2utm)
    # Make and save a colored 'lakes' overlay.
    ice_elevation = 1000.0 # Hardcoded that lakes above 1000m are frozen.
    save_lakes_overlay_png(segments, elevation, cell_iter, ice_elevation, fofo)
    true
end

function save_lakes_overlay_png(segments, elevations, cell_iter, ice_elevation, folder)
    # Create a black-and-white image, where any segments with any trace of a lake is a lake.
    # This step is unnecessary legacy from an older function structure, refactoring incomplete.
    # There's no need to make this matrix now.
    lm_bool = map(labels_map(segments)) do i
        islake = segment_mean(segments, i) > 0.0f0
        Gray{Bool}(islake)
    end
    # Hardcoded lake colors
    water_color = RGBA{N0f8}(0.521, 0.633, 0.764, 1.0)
    ice_color = RGBA{N0f8}(0.8, 0.8, 0.8, 1.0)
    transparent = RGBA{N0f8}(0.0, 0.0, 0.0, 0.0)
    # Create the colourful, transparent image (in full resolution, disregarding cell_to_utm_factor)    
    img = map(zip(lm_bool, elevations[cell_iter])) do (is_lake, elevation)
        if is_lake == Gray{Bool}(true)
            if elevation > ice_elevation
                ice_color
            else
                water_color
            end
        else
            transparent
        end
    end
    # Feedback
    display_if_vscode(img)
    # Save
    ffna = joinpath(folder, WATER_FNAM)
    @debug "    Saving $ffna"
    save_png_with_phys(ffna, img)
    # Feedback for testing purposes
    img
end

"""
    lake_segment(elevations, cell_iter, horizontal_distance, lake_steepness_max)
    ---> Segmented image
"""
function lake_segment(elevations, cell_iter, horizontal_distance, lake_steepness_max)
    # Hardcoded parmeter for identifying lake regions from steepness using Felzenszwalb regions
    lake_area_min = 900 # m²
    lake_pixels_min = Int(round(lake_area_min / horizontal_distance^2))
    k = 3.8
    @debug "    For water id, finding local steepness over 2 m lenghts"
    steepness = steepness_decirad_capped(elevations, cell_iter, horizontal_distance)
    # Boolean result matrix
    @debug "    For water id, classifying areas of low steepness"
    _lake_segment(steepness, k, lake_steepness_max, lake_pixels_min)
end


"""
    steepness_decirad_capped(elevations, cell_iter)
    ---> Matrix{Float32}

`elevation` is a matrix with grid distance 1 unit, same unit as elevation.

Returns local stepness around each pixel, a matrix of scalars in units of decirad:

    10π / 2 decirad == 90°
    10π / 180 decirad == 1°
    1 decirad == 180 / 10π ° == 5.729577951308232°
"""
function steepness_decirad_capped(elevation, cell_iter, cell2utm)
    # Output image size
    ny, nx = size(cell_iter)
    # Source indices
    indices = (1:cell2utm:(ny * cell2utm), 1:cell2utm:(nx * cell2utm))
    # Apply
    mapwindow(local_steepness_capped, elevation, (3, 3); indices)
end

@inbounds function local_steepness_capped(M::Matrix)
    # This function will operate on the original data resolution,
    # so as to be independent of output resolution, cell to utm factor.
    # This will find the local steepness around each output pixel.
    @assert size(M) == (3, 3)
    # If rows in M correspond to south -> north
    # and cols in M correspond to west -> east
    # _ n _
    # w z e  
    # _ s _
    _, w, _, n, z, s, _, e, _ = M
    gr1 = (e - w) / 2
    gr2 = (n - s) / 2
    # We're capping steepness at
    cap_deg =  1.5f0
    # For low steepness, atan(x) == x
    10 * min(hypot(gr1, gr2), cap_deg * π / 180)
end



"""
    _lake_segment(steepness, k, lake_steepness_max, lake_pixels_min)
    ---> Segmented image
"""
function _lake_segment(steepness, k, lake_steepness_max, lake_pixels_min)
    # Divide the steepness matrix into segments, based on proximity and steepness difference.
    # Note that we don't have or use a minimum cell count argument yet.
    # Note: A size check might be good here because this can take > 10 minutes.
    @debug "    First steepness segmentation"
    steep_segments = felzenszwalb(steepness, k)
    is_flat(i) = segment_mean(steep_segments, i) < lake_steepness_max
    is_large(i) = segment_pixel_count(steep_segments, i) >= lake_pixels_min
    # Create a black-and-white image, where we discard small and too steep segments.
    # We do not create a bool matrix, because we want to re-segment the result afterwards.
    @debug "    Grouping segments"
    islake_matrix = map(labels_map(steep_segments)) do i # i is a label, representing a set of pixels.
        Gray{N0f8}(is_large(i) && is_flat(i))
    end
    # A common artifact is lines across a lake, possibly from power lines.
    # Let's grow the lakes, then shrink, with 8-connectivity to add
    # lake shores.
    @debug "    Dilate -> erode lakes"
    dilate!(islake_matrix, copy(islake_matrix))
    erode!(islake_matrix, copy(islake_matrix))
    # We are confident that positives are true positives. Still, we have lots of fake negatives
    # inside of the lake regions. They would appear as islands that are really just noise, waves,
    # or recalibration. Most of them coindicentally fall below lake_pixels_min:
    @debug "    Second segmentation"
    lake_segments = felzenszwalb(islake_matrix, k, lake_pixels_min)
    #=
    # Check that the first segment is NOT a lake.
    i_1 = findfirst(labels_map(lake_segments)) do l 
        l == 1
    end
    is_1_lake = segment_mean(lake_segments, i_1) < lake_steepness_max && isegment_pixel_count(steep_segments, i) >= lake_pixels_min
    throw("TODO devise a working test")
    =#
end


function sum_of_cartesian_indices_per_label(segments)
    R = CartesianIndices(axes(labels_map(segments)))
    dic = Dict{Int64, CartesianIndex}()
    for ci in R
        label = labels_map(segments)[ci]
        if haskey(dic, label)
            dic[label] += ci
        else
            push!(dic, label => ci)
        end
    end
    dic
end

function centre_indices(segments)
    R = CartesianIndices(axes(labels_map(segments)))
    diccum = sum_of_cartesian_indices_per_label(segments)
    dicmean = Dict{Int64, CartesianIndex}()
    for label in keys(diccum)
        cum = diccum[label]
        count = segment_pixel_count(segments, label)
        mean = CartesianIndex{2}((cum.I .÷ count)) # We stick to integer indices here.
        if mean ∉ R
            @show mean cum count
            throw("what?")
        end 
        push!(dicmean, label => mean)
    end
    dicmean
end

function write_lake_to_csv(ffnam_csv_lakes, lakes::SegmentedImage, elevation, cell2utm)
    headers = ["Elevation_m", "Area_m²", "Sheet_index"]
    # 
    vz = Int[]
    vA = String[]
    vij = String[]
    # Indices of the centre of area for each lake
    dic = centre_indices(lakes)
    for (label, center_i) in dic
        if ! iszero(segment_mean(lakes, label))
            z = Int(round(elevation[cell2utm * center_i]))
            A = segment_pixel_count(lakes, label) * cell2utm^2
            push!(vz, z)
            push!(vA, string(A))
            push!(vij, string(center_i.I))
        end
    end
    # Order by surface area
    order = sortperm(vz; rev = true)
    vectors = [string.(vz)[order], vA[order], vij[order]]
    widths = [20, 20, 20]
    write_vectors_to_csv(ffnam_csv_lakes, headers, vectors, widths)
end