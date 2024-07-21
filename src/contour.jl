# Step in pipeline.
# Creates elevation contour lines overlay from elevation data.
# # Output is an image file per sheet, for manual touch-up.

"""
    contour_lines_overlay(sb::SheetBuilder)
    contour_lines_overlay(fofo, cell_iter, cell_to_utm_factor)

"""
function contour_lines_overlay(sb::SheetBuilder)
    # Line thickness definitions are not part of the SheetBuilderMatrix.
    # Instead, we read those directly from the .ini file. (This may
    # possibly stand in the way of parallelisation for large matrices.)
    thick_20 = get_config_value("Line thicknesses", "Elevation contour 20m", Int; nothing_if_not_found = false)
    thick_100 = get_config_value("Line thicknesses", "Elevation contour 100m", Int; nothing_if_not_found = false)
    thick_1000 = get_config_value("Line thicknesses", "Elevation contour 1000m", Int; nothing_if_not_found = false)
    contour_lines_overlay(full_folder_path(sb), sb.cell_iter, cell_to_utm_factor(sb), thick_20, thick_100, thick_1000)
end
function contour_lines_overlay(fofo, cell_iter, cell2utm, thick_20, thick_100, thick_1000)
    if isfile(joinpath(fofo, COMPOSITE_FNAM))
        @debug "    $COMPOSITE_FNAM in $fofo already exists. Exiting `contour_lines_overlay`"
        return true
    end
    if isfile(joinpath(fofo, CONTOUR_FNAM))
        @debug "    $CONTOUR_FNAM in $fofo already exists. Exiting `contour_lines_overlay`"
        return true
    end
    res = _elev_contours(fofo, cell_iter, cell2utm, thick_20, thick_100, thick_1000)
    # Feedback
    display_if_vscode(res)
    # Save
    ffna = joinpath(fofo, CONTOUR_FNAM)
    @debug "    Saving $ffna"
    save_png_with_phys(ffna, res)
    true
end
function _elev_contours(fofo, cell_iter, cell2utm, thick_20, thick_100, thick_1000)
    #################
    # Allocate arrays
    #################
    # Read elevation
    g = readclose(joinpath(fofo, CONSOLIDATED_FNAM))
    # Prepare iterators and buffers
    ny, nx = size(cell_iter)
    # A better way would likely be to transpose the source data.
    # However, this works, too. 
    source_indices = (1:cell2utm:(nx  * cell2utm), 1:cell2utm:(ny * cell2utm))
    si = CartesianIndices(source_indices)
    # Store elevation in an 'image buffer'. We'll add more below.
    # We store info in an 'image' in order to use the fast image filtering functions later.
    zs = GrayA{Float32}.(g.A[:, :, 1], Array{Float32}(undef, size(g.A[:, :, 1])...))
    # Store nearby steepness in the alpha channel of the 'image buffer'.
    channelview(zs)[2, :, :] .= terrain_steepness(zs)
    # Pre-allocate boolean buffer. This has the same size as the output image
    bbuf = Array{Gray{Bool}}(undef, size(si)...)
    # Pre-allocate boolean result image
    res = zeros(Gray{Bool}, size(si)...)
    # With all buffers ready, call the inner function
    _elev_contours!(res, source_indices, zs, bbuf, [thick_20, thick_100, thick_1000], [20, 100, 1000])
    # Go from black-and-white to defined colours. Flip axes to image-like.
    transpose(map(res) do pix
        pix == true && return RGBA{N0f8}(0.173, 0.192, 0.255, 1)
        RGBA{N0f8}(0., 0, 0, 0)
    end)
end
function _elev_contours!(res, source_indices, zs, bbuf, thicknesses, elevation_spacings)
    for (t, Δz) in zip(thicknesses, elevation_spacings)
        # Overwrite bbuf with pixels on contours.
        mapwindow!(func_elev_contour(Δz), bbuf, zs, (1, 1); indices = source_indices)
        # We have tried to reduce contours on bumps. Now remove most of the rest:
        # Our very clever thinning, despeckling and thickening:
        bbuf .= strokeify(bbuf, t)
        # Treat 'false' as transparent. Overlay bbuf on res and write in-place to res.
        map!(BlendLighten, res, res, bbuf)
    end
    res
end

function terrain_steepness(zs)
    # Prepare a wide low-pass filter. Window size:
    w = 49
    # Coefficients, including a window.
    c = Float32.(fir_lp_coefficients(w))
    # Elevations smoothed by lp-filter
    zf = imfilter(channelview(zs)[1, :, :], (c, transpose(c)), FIRTiled());
    # Find nearby steepness, i.e. not the steepness of the
    # local stone, tree or shed.
    hypot.(imgradients(zf, KernelFactors.ando5, "replicate")...)
end


"""
    strokeify(bw, thickness)

Applies a series of filters that change a collection of adjacent pixels to
pen strokes with the given thickness in pixels. Less exact for diagonal strokes.
"""
function strokeify(bw, thickness)
    thickness == 0 && return bw
    if (thickness - 1 ) % 2 !== 0
        throw(ArgumentError("thickness = $thickness ∉ [0, 1, 3, ...] "))
    end
    # Despeckle
    remove_isolated_pixels!(bw)
    # Reduce all the 'true' regions to one pixel width.
    stroked = reinterpret(Gray{Bool}, thinning(reinterpret(Bool, bw), algo = GuoAlgo()))
    # Remove pixel islands, 5 or smaller
    remove_small_islands!(stroked, 5)
    # Make all the thin regions as thick as specified
    thickness == 1 && return stroked
    dilate(stroked; r = thickness ÷ 2)
end

func_elev_contour(vert_dista) = func_elev_contour(Float32(vert_dista))
function func_elev_contour(vert_dista::Float32)
    # Hardcoded, but optimum likely depends on terrain types and measurement.
    # A value of 0.0 minimises isolated pixels; 0.6 emphasises small features.
    noise_error = 0.0f0
    # This works fine with value 1.0. Increasing it yields more continuous lines in steep terrain.
    steepness_allowance = 2.0
    (M) -> begin
        pix = first(M)
        z = gray(pix)
        steepness = alpha(pix)
        # No elevation contour line at sea level
        z < vert_dista && return Gray{Bool}(false)
        # This cell is on a contour if the sampled elevation level
        # is within an interval. The interval size depends on
        # steepness in the proximity, but not on local steepness, which is
        # partly due to large stones, trees, houses.
        # In this way, a cell on a house, stone, or a tree is less likely
        # to be deemed on an elevation contour.
        if z % vert_dista < (steepness_allowance * steepness + noise_error)
            return Gray{Bool}(true)
        end
        Gray{Bool}(false)
    end

end
function black_white_elevation_contour!(bw, z, elevation, thick_1000)
    bw
end

function remove_isolated_pixels!(img::Array{Gray{Bool},2})
    kernel = centered(Int[1 1 1; 1 0 1; 1 1 1])
    # Make a similar matrix{Float32} counting neighbors of a pixel
    neighbors = imfilter(img, kernel, Fill(Gray{Bool}(false)))
    # Very clever: ´x .&= y´ is a synonym for ´x .= x .& y´
    #               ´&´ is bitwise ´and´.
    # So, for every pixel:
    #      1) If a pixel is not set, set it to false without considering neigbors.
    #      2) If a pixel is set (x = true), and it does not have neighbors (y is false),
    #         unset it.
    #      3) If a pixel is set, and it does have neighbors, set it.
    reinterpret(Bool, img) .&= (neighbors .> 0)
    img
end

# TODO use this from 'water'
function remove_small_islands!(img, max_pixels)
    segments = felzenszwalb(reinterpret(Bool, img), 1.0, 2)
    labels = labels_map(segments)
    # A dictionary of pixels in each segment no.
    segment_counts = segment_pixel_count(segments)
    # Create a boolean array for valid segments (larger than max_pixels and nonblank)
    valid_segments = Dict(i => (segment_counts[i] > max_pixels && segment_mean(segments, i) == 1) for i in keys(segment_counts))
    # Update the img array in-place
    for i in eachindex(img)
        img[i] = Gray{Bool}(valid_segments[labels[i]])
    end
    img
end
