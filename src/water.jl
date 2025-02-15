# Step in pipeline.
# Identifying water surfaces is a time consuming operation.
# Results look good, but may benefit from manual touch-up of the output file.
# This might probably be speed up a good deal.
"""
    water_overlay(sb::SheetBuilder)
    water_overlay(fofo, cell_iter, cell2utm, f_I_to_utm)
    ---> Bool

Creates water overlay from elevation data.
Output is an RGBA .png image file. Fully transparent outside water surfaces, for manual touch-up.
"""
function water_overlay(sb::SheetBuilder)
    # Go ahead
    water_overlay(full_folder_path(sb), sb.cell_iter, cell_to_utm_factor(sb), sb.f_I_to_utm)
end
function water_overlay(fofo, cell_iter, cell2utm, f_I_to_utm)
    # Early exits
    if ! isfile(joinpath(fofo, CONSOLIDATED_FNAM))
        @debug "    $CONSOLIDATED_FNAM in $fofo\n           does not exist. Exiting `water_overlay`"
        return false
    end
    ffnam_csv_lakes = splitext(joinpath(fofo, WATER_FNAM))[1] * ".csv"
    if isfile(joinpath(fofo, WATER_FNAM)) && isfile(ffnam_csv_lakes) 
        @debug "    $WATER_FNAM and $(splitpath(ffnam_csv_lakes)[end]) in $fofo\n           already exists. Exiting `water_overlay`"
        return true
    end
    # Harvest defaults. 
    sea_level = get_config_value("Water detection", "Sea level max [m]", Float32)
    area_lim = get_config_value("Water detection", "Minimum area [m²]", Int64)
    slopes_local_lim = get_config_value("Water detection", "Local slopes [rad] limit", Float32)
    r_slopes = get_config_value("Water detection", "Radius [m], standard dev. for slopes blurring", Int64)
    Δzlim = get_config_value("Water detection", "Elevation range [m] limit for a water body", Float32)
    diagonal_length_lim = get_config_value("Water detection", "Bounding box diagonal length [m] limit", Float64)
    r_art = get_config_value("Water detection", "Radius [m] for artifact criterion blurring", Int64)
    artifact_lim = get_config_value("Water detection", "Artifact [m⁻²] limit", Float64)
    flux_lim = get_config_value("Water detection", "Flux [m] limit", Float64)
    water_overlay(fofo, cell_iter, cell2utm, f_I_to_utm,
        sea_level, area_lim, slopes_local_lim, 
        r_slopes, Δzlim, diagonal_length_lim, 
        r_art, artifact_lim, flux_lim)
    true
end
function water_overlay(fofo, cell_iter, cell2utm, f_I_to_utm,
    sea_level, area_lim, slopes_local_lim, 
    r_slopes, Δzlim, diagonal_length_lim, 
    r_art, artifact_lim, flux_lim)
    # Find segments with lakes, and without. Whether a segment is water
    # can be determined by: `i -> segment_mean(segments, i) > 0.5`
    segments = water_segments(fofo, cell2utm,
        sea_level, area_lim, slopes_local_lim, 
        r_slopes, Δzlim, diagonal_length_lim, 
        r_art, artifact_lim, flux_lim)
    @debug "    Water overlay, output .csv and .png"
    # Identify water segment label numbers (the remaining labels are terrain)
    # This includes sea, of course
    iswater = i -> segment_mean(segments, i) > 0.5 
    labels_water = filter(iswater, segment_labels(segments))
    # Indices of all segments in output image
    index_dic = mean_index_dictionary(segments)
    # Output resolution elevations 
    z = elevation_at_output(fofo, cell_iter, cell2utm)
    # Segment elevations dictionary
    z_dic = Dict([label => z[I] for (label, I ) in index_dic])
    # Indices of water segments
    vI = [index_dic[label] for label in labels_water]
    # Segment area
    va = [cell2utm^2 * segment_pixel_count(segments, label) for label in labels_water]
    # Geographical utm positions of water surfaces
    vutm = f_I_to_utm.(vI)
    # Vector of all water segment elevations
    vz = [Int64(round(z_dic[label])) for label in labels_water]
    # Save data in a .csv file
    write_water_surfaces_to_csv(fofo, vz, va, vutm, getfield.(vI, :I))
    # Make and save a colored 'lakes' overlay.
    save_water_overlay_png(fofo, segments, z_dic)
end

function write_water_surfaces_to_csv(fofo, vz, va, vutm, vItup)
    ffnam_csv_lakes = splitext(joinpath(fofo, WATER_FNAM))[1] * ".csv"
    @assert length(vz) == length(va) == length(vutm) == length(vItup)
    nt = (;Elevation_m = vz, Area_m² = va, Utm = vutm, Cell_index = vItup)
    write_named_tuple_to_csv(ffnam_csv_lakes, nt)
    true
end    


function save_water_overlay_png(fofo, segments, z_dic)
    iswater = i -> segment_mean(segments, i) > 0.5
    # Hardcoded lake colors
    water_color = RGBA{N0f8}(0.521, 0.633, 0.764, 1.0)
    ice_color = RGBA{N0f8}(0.8, 0.8, 0.8, 1.0)
    transparent = RGBA{N0f8}(0.0, 0.0, 0.0, 0.0)
    # Colourful image
    iswater = label -> segment_mean(segments, label) > 0.5
    isicewater = label -> z_dic[label] > 1000
    img = map(labels_map(segments)) do label
        if iswater(label)
            if isicewater(label)
                ice_color
            else
                water_color
            end
        else
            transparent
        end
    end
    # Feedback
    # display_if_vscode(img)
    # Save
    ffna = joinpath(fofo, WATER_FNAM)
    @debug "    Saving $ffna"
    save_png_with_phys(ffna, img)
    # Feedback for testing purposes
    img
end


function water_segments(fofo, cell2utm,
        sea_level, area_lim, slopes_local_lim, 
        r_slopes, Δzlim, diagonal_length_lim, 
        r_art, artifact_lim, flux_lim)
    z = elevation_full_res(fofo)
    ∇²z = divergence_of_gradients(z)
    @debug "    Water overlay, find candidate segments"
    candidate_segments = water_candidate_segments(cell2utm, z, ∇²z,
        sea_level, area_lim, slopes_local_lim, 
        r_slopes, Δzlim, diagonal_length_lim, 
        r_art, artifact_lim)
    # Feedback
    # display_if_vscode(candidate_segments; randomcolor = true)
    @debug "    Water overlay, prune by area flux"
    # Basis on which to eliminate a segment. We're not finding actual flux here, but 
    # a downscaled version of it. See `flux_lim_scaled`.
    dic = flux_of_segments_dictionary(downsample(∇²z, cell2utm), candidate_segments)
    # Criterion for removing a candidate segment. Scaling this is interesting,
    # since flux is defined along the perimeter of a segment, which scales linearly,
    # but we use Green's theorem on a number of pixels with divergence of gradients.
    # So we scale the criterion with area instead of linearly.
    flux_lim_scaled = flux_lim / cell2utm^2 
    # Candidate segment means are not just 0 or 1 because of regions merging in the previous step. 
    f_to_be_removed = i -> segment_mean(candidate_segments, i) > 0.5 && dic[i] <= flux_lim_scaled
    # In case a segment to be removed has neighbors, it will be merged with the one minimizing this function:
    f_diff = (rem_label, neigh_label) -> segment_mean(candidate_segments, neigh_label)
    # Merge segments with too low flux (i.e. that are taller than their surroundings)
    segments = prune_segments(candidate_segments, f_to_be_removed, f_diff)
    # Feedback
    display_if_vscode(segments; randomcolor = true)
    segments
end
function water_candidate_segments(cell2utm, z, ∇²z,
    sea_level, area_lim, slopes_local_lim, 
    r_slopes, Δzlim, diagonal_length_lim, 
    r_art, artifact_lim)
    #
    # We have three components, one of which should be true for a cell to be water.
    # is_flat_large_water has been checked for size, but not the two other ones.
    # (E.g. an artifact is allowed to be smaller, as it typically covers part of a water surface)
    is_water_candidate  = is_local_and_bbox_diagonal_slope_in_range(cell2utm, z, slopes_local_lim, r_slopes, Δzlim, diagonal_length_lim) .| 
        is_artifact(cell2utm, ∇²z, r_art, artifact_lim) .| 
        is_low_water(cell2utm, z, sea_level)
    # Feedback
    # display_if_vscode(is_water_candidate)

    # Segment the boolean matrix, while applying size limit.  
    # 'k' is not sensitive when there is only two values like here.
    # The result, since some segments are joined from '1' and '0', does not 
    # exclusively contain  mean values 0 or 1.
    @debug "    Water overlay, join candidates"
    pixel_lim = Int(round(area_lim / cell2utm^2))
    segments = felzenszwalb(is_water_candidate, 6, pixel_lim)
    # Feedback
    # display_if_vscode(segments; randomcolor = true)
    segments
end

function is_local_and_bbox_diagonal_slope_in_range(cell2utm, z, slopes_local_lim, r_slopes, Δzlim, diagonal_length_lim)
    # We fill sections with unacceptable local steepness with this:   
    no_water_value = Float32(2 * slopes_local_lim)
    # Segmented image, downsampled by cell2utm
    segments = locally_acceptable_steepness_segments(cell2utm, z, no_water_value, slopes_local_lim, r_slopes)
    # Feedback
    # display_if_vscode(segments; randomcolor=true)
    is_local_flat(i) = segment_mean(segments, i) <= no_water_value
    # Function `is_local_flat` identifies some reasonably flat segments, probably including agricultural fields with minor slopes. 
    # We will put additional demands on the vertical variation within each segment, i.e. lake. For small lakes, we allow propotionally less variation.
    # The actual sea data won't necessarily meet all of these criteria, because the sea may have elevation differences due to tides.
    # We will make exceptions for the sea elsewhere. The same goes for artifacts which are part of lakes.
    #
    # Dictionaries for max and min elevation within each segment
    min_z_dic = segment_dictionary(downsample(z, cell2utm), segments, min)
    max_z_dic = segment_dictionary(downsample(z, cell2utm), segments, max)
    # Dictionary for each segment's bounding box diagonal (kind of like the enclosing diameter, but dependant  on coordinate system orientation)
    bb_diag_len_dic = bbox_diagonal_length_dictionary(segments)
    function is_bbox_diagonal_slope_in_range(i)
        maxz = max_z_dic[i]
        minz = min_z_dic[i]
        diagonal_length = min(bb_diag_len_dic[i], diagonal_length_lim / cell2utm)
        maxz - minz <= Δzlim * diagonal_length * cell2utm / diagonal_length_lim
    end
    #=
    # DEBUG
    dicind = mean_index_dictionary(segments)
    for i in segment_labels(segments)
        print(rpad("$i", 5))
        print(rpad("Centre $(dicind[i].I)", 30))
        print(rpad("Pixel count $(segment_pixel_count(segments, i))", 25))
        print(rpad("z $(round(max_z_dic[i])) m", 20))
        print(rpad("Δz $(round(max_z_dic[i] - min_z_dic[i], digits=2))", 10))
        diagonal_length = min(bb_diag_len_dic[i], diagonal_length_lim)
        print(rpad("diagonal_length for slope: $(round(diagonal_length))", 40))
        println(rpad("slope in range: $(is_bbox_diagonal_slope_in_range(i))", 20))
    end
    # /DEBUG
    =#
    BitMatrix(map(labels_map(segments)) do i # i is the label of the current pixel in segments
        is_local_flat(i) && is_bbox_diagonal_slope_in_range(i)
    end)
end
is_low_water(cell2utm, z, sea_level) = downsample(z, cell2utm) .<= sea_level
function is_artifact(cell2utm, ∇²z, r_art, artifact_lim)
    wz = weighted_∇⁴z(∇²z)
    # wz is an indicator of artifacts presence at a cell,
    # but we need to blur it to find contiguous areas of artifacts
    # and drop the mostly isolated false identifications.
    downsample(imfilter(wz, Kernel.gaussian(r_art)), cell2utm) .> artifact_lim
end

function weighted_∇⁴z(∇²z; x₀::Float64 = 10^-6, x₁::Float64 = 10^(-3.5))
    # ∇²z is the divergence of the gradient of scalar field z. 
    # This divergence is of course a scalar field.
    # We're looking at how the divergence (i.e. curvature) changes with space,
    # in other words the bi-laplacian ∇²(∇²z).
    ∇⁴z = .+(jacobian_components(∇²z)...)
    # We have found that in areas with artifacts, values of ∇⁴z 
    # is mostly in the range (x₀, x₁). However, 
    # the values overlap with some 'naturally occuring' values.
    # 
    w = func_tanh_window(x₀, x₁)
    # So these values will be close to 1 where artifacts are present:
    w.(abs.(∇⁴z))
end

"""
    func_tanh_window(x₀, x₁)

If x₀ < x₁, returns a window function with maximum 1.0, sloping off symmetrically on a log10 axis.
Also see the one-sided version, `func_tanh_limit`.

Example
```
julia> f = func_tanh_window(1e-1, 1)
#16 (generic function with 1 method)

julia> vx = log10_rng(1e-5, 1e4, 10);

julia> hcat(vx, f.(vx))
10×2 Matrix{Float64}:
     1.0e-5  0.000627443
     0.0001  0.00462496
     0.001   0.0335707
     0.01    0.219028
     0.1     0.824027
     1.0     0.824027
    10.0     0.219028
   100.0     0.0335707
  1000.0     0.00462496
 10000.0     0.000627443

julia> f(1e-1)
0.8240271368319426

julia> f(1)
0.8240271368319426

julia> f(√(1e-1 * 1))
1.0
```

"""
function func_tanh_window(x₀::Float64, x₁::Float64)
    @assert x₀ < x₁
    # The maximum is, not so obviously, at √(x₀x₁)
    # We're scaling so as to have maximum 1.0.
    let ymax = (1 / 4 ) * (1 + tanh(log10(√(x₀ * x₁) / x₀))) * (1 - tanh(log10(√(x₀ * x₁) / x₁)))
       x -> (1 / 4 ) * (1 + tanh(log10(x / x₀))) * (1 - tanh(log10(x / x₁))) / ymax
    end
end

function locally_acceptable_steepness_segments(cell2utm, z, no_water_value, slopes_local_lim, r_slopes)
    # Where steepness is above the slopes_local_lim, we set
    # the value to 'no_water_value', for easy segmentation below.
    st = local_steepness_where_acceptable(z, no_water_value, slopes_local_lim, r_slopes)
    k = 6
    # We set a minimum number of pixels here for speed, but
    # larger segments are filtered out elsewhere, depending on the configureable "Minimum area [m²]"
    preliminary_minimum_pixels_count = Int(round(200 / cell2utm^2))
    felzenszwalb(downsample(st, cell2utm), k, preliminary_minimum_pixels_count)
end

function local_steepness_where_acceptable(z, no_water_value::Float32, slopes_local_lim::Float32 = 0.029f0, r_slopes = 2)
    ssteep = Float32.(smoothed_steepness(z, r_slopes))
    map(ssteep) do st 
        st .<= slopes_local_lim ? st : no_water_value
    end
end
function smoothed_steepness(z, r_slopes)
    # We're smoothing before calculating steepness, because
    # the derivatives have sign (whereas steepness does not have sign)
    hypot.(blurred_gaussian(z, r_slopes)...)
end

function blurred_gaussian(z, r_slopes)
    z1, z2 = jacobian_components(z)
    imfilter(z1, Kernel.gaussian(r_slopes)), imfilter(z2, Kernel.gaussian(r_slopes))
end
function flux_of_segments_dictionary(∇²z, segments)
    # Flux for all segments, valid or not
    flux_dict = segment_dictionary(∇²z, segments, +)
    # At the border of the elevation matrix, the divergence simply reflects 
    # the general border conditions we used for the calculation.
    # Hence, we invalidate the flux of a segment if the area is bordering the matrix,
    # by setting the flux of that segment to zero. 
    R = CartesianIndices(axes(labels_map(segments)))
    segment_bbox_dict = bbox_internal_indices_dictionary(segments)
    for label in segment_labels(segments)
        bb = segment_bbox_dict[label]
        tl = bb[1] + CartesianIndex(-1, -1)
        br =  bb[end] + CartesianIndex(1, 1)
        if tl ∉ R || br ∉ R
            # Stepping one index away from at least one corner brought us outside the matrix.
            # We don't trust this flux calculation.
            flux_dict[label] = zero(flux_dict[label])
        end
    end
    # Flux for valid segments, zero for segments at border
    flux_dict
end