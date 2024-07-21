# Step in pipeline.
# Identifying water surfaces is a time consuming operation.
# Results look good, but may benefit from manual touch-up of the output file.
# This might probably be speed up a good deal.
"""
    water_overlay(sb::SheetBuilder)
    water_overlay(fofo)

Creates water overlay from elevation data.
Output is image files with full transparency outside water surfaces, for manual touch-up.
"""
function water_overlay(sb::SheetBuilder)
    lake_steepness_max = get_config_value("Water", "Lake steepness max", Float32; nothing_if_not_found = false)
    # Go ahead
    water_overlay(full_folder_path(sb), sb.cell_iter, cell_to_utm_factor(sb), lake_steepness_max)
end
function water_overlay(fofo, cell_iter, cell2utm, lake_steepness_max)
    if isfile(joinpath(fofo, COMPOSITE_FNAM))
        @debug "    $COMPOSITE_FNAM in $fofo already exists. Exiting `water_overlay`"
        return true
    end
    if ! isfile(joinpath(fofo, CONSOLIDATED_FNAM))
        @debug "    $CONSOLIDATED_FNAM in $fofo does not exist. Exiting `water_overlay`"
        return false
    end
    if isfile(joinpath(fofo, WATER_FNAM))
        @debug "    $WATER_FNAM in $fofo already exists. Exiting `water_overlay`"
        return true
    end
    ffna = joinpath(fofo, CONSOLIDATED_FNAM)
    elevations = let
        z = readclose(ffna)
        z.A[:, :, 1]
    end
    lm_bool = is_water_surface(elevations, cell_iter, cell2utm, lake_steepness_max)
    ice_elevation = 1000.0 # Hardcoded that lakes above 1000m are frozen.
    save_lakes_overlay_png(lm_bool, elevations, ice_elevation, fofo)
    true
end

function save_lakes_overlay_png(lm_bool, elevations, ice_elevation, folder)
    # Hardcoded lake colors
    water_color = RGBA{N0f8}(0.521, 0.633, 0.764, 1.0)
    ice_color = RGBA{N0f8}(0.8, 0.8, 0.8, 1.0)
    transparent = RGBA{N0f8}(0.0, 0.0, 0.0, 0.0)
    # Create the colourful, transparent image (in full resolution, disregarding cell_to_utm_factor)
    # TODO Speed up, see contour.jl
    img = transpose(map(zip(lm_bool, elevations)) do (is_lake, elevation)
        if is_lake == Gray{Bool}(true)
            if elevation > ice_elevation
                ice_color
            else
                water_color
            end
        else
            transparent
        end
    end)
    # Feedback
    display_if_vscode(img)
    # Save
    ffna = joinpath(folder, WATER_FNAM)
    @debug "    Saving $ffna"
    save_png_with_phys(ffna, img)
    # Feedback for testing purposes
    img
end


function is_water_surface(elevations, cell_iter, horizontal_distance, lake_steepness_max)
    # Hardcoded parmeter for identifying lake regions from steepness using Felzenszwalb regions
    lake_area_min = 900 # m²
    lake_pixels_min = Int(round(lake_area_min / horizontal_distance^2))
    k = 3.8
    @debug "    Calculating steepness"
    steep_matrix = steepness_decirad_capped(elevations, cell_iter, horizontal_distance)
    # Boolean result matrix
    lake_matrix(steep_matrix, k, lake_steepness_max, lake_pixels_min)
end


"""
    steepness_decirad_capped(elevations, cell_iter, horizontal_distance; cap_deg = 1.5)
    ---> Matrix{Float64}

`elevations` is a matrix in the same unit as scalar `horizontal_distance`.

Returns local stepness, a matrix of scalars in values of decirad:

    10π / 2 decirad == 90°
    10π / 180 decirad == 1°
    1 decirad == 180 / 10π ° == 5.729577951308232°

For other purposes than finding water surfaces, cap_deg = 86 ° is a practical upper limit.
"""
function steepness_decirad_capped(elevations, cell_iter, horizontal_distance; cap_deg =  1.5)
    # WORK IN PROGRESS.
    # Try the same approach as in topo_relief.
    # Transpose the data before calling this.
    # Don't allocate g1, g1, use map or mapwindow instead.
    z_norm = elevations ./ horizontal_distance
    # The next line may meet memory restrictions. It calculates gradients at every sample point.
    # If we sampled fewer points, the parameters for dealing with elevation data noise 
    # would need adjustments.
    # One workaround may be to split areas further, and patch together 'WATER_FNAM' manually.
    g1, g2 = imgradients(z_norm, KernelFactors.prewitt)
    # We limit inclination to 86 degree, as a compromise. Many above that value are artifacts.
    # TODO speed up, see contours.jl
    map(zip(g1, g2)) do (gr1, gr2)
        min(10 * atan(hypot(gr1, gr2)), cap_deg * 10 * π / 180)
    end
end



function lake_matrix(steep_matrix, k, lake_steepness_max, lake_pixels_min)
    # Divide the steepness matrix into segments, based on proximity and steepness difference.
    # Note that we don't have or use a minimum cell count argument here. 
    # Note: A size check might be good here because this can take > 10 minutes.
    @debug "    First steepness segmentation"
    steep_segments = felzenszwalb(steep_matrix, k)
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
    cleanup_segments = felzenszwalb(islake_matrix, k, lake_pixels_min)
    # Create a black-and-white image, where any segments with any trace of a lake is a lake.
    @debug "    Cleaning fake islands"
    map(labels_map(cleanup_segments)) do i
        islake = segment_mean(cleanup_segments, i) > 0.0f0
        Gray{Bool}(islake)
    end
end

