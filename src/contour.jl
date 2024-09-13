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
    thick_20 = get_config_value("Elevation contour lines", "Thickness 20m", Int)
    thick_100 = get_config_value("Elevation contour lines", "Thickness 100m", Int)
    thick_1000 = get_config_value("Elevation contour lines", "Thickness 1000m", Int)
    vthick = [thick_20, thick_100, thick_1000]
    vdist = [20, 100, 1000] # Contour spacings.
    minlen = get_config_value("Elevation contour lines", "Minimum length", Int)
    contour_lines_overlay(full_folder_path(sb), sb.cell_iter, cell_to_utm_factor(sb), minlen, vthick, vdist)
end
function contour_lines_overlay(fofo, cell_iter, cell2utm, minlen, vthick, vdist)
    if isfile(joinpath(fofo, CONTOUR_FNAM))
        @debug "    $CONTOUR_FNAM in $fofo \n           already exists. Exiting `contour_lines_overlay`"
        return true
    end
    if ! isfile(joinpath(fofo, CONSOLIDATED_FNAM))
        @debug "    $CONSOLIDATED_FNAM in $fofo\n           does not exist. Exiting `contour_lines_overlay`"
        return false
    end
    res = _elev_contours(fofo, cell_iter, cell2utm, minlen, vthick, vdist)
    # Feedback
    #display_if_vscode(res)
    # Save
    ffna = joinpath(fofo, CONTOUR_FNAM)
    @debug "    Saving $ffna"
    save_png_with_phys(ffna, res)
    true
end
function _elev_contours(fofo, cell_iter, cell2utm, minlen, vthick, vdist)
    #################
    # Allocate arrays
    #################
    # Read elevation
    g = readclose(joinpath(fofo, CONSOLIDATED_FNAM))
    #
    # Prepare iterators and buffers, reduce source resolution
    #
    ny, nx = size(cell_iter)
    source_indices = (1:cell2utm:(ny  * cell2utm), 1:cell2utm:(nx * cell2utm))
    si = CartesianIndices(source_indices)
    # After this resolution reduction we can use cell_iter for both source and dest.
    elevation = transpose(g.A[:, :, 1])[si]
    y´, x´ = terrain_gradient(elevation)
    # We elevation and gradient info in an 'image' in order to use the fast image filtering functions later.
    zs = RGB{Float32}.(elevation, y´, x´)
    # Pre-allocate boolean buffer (we make one countour distance at a time)
    bbuf = Array{Gray{Bool}}(undef, size(si)...)
    # Pre-allocate boolean result image
    res = zeros(Gray{Bool}, size(si)...)
    # With all buffers ready, call the inner function
    _elev_contours!(res, zs, bbuf, minlen, vthick, vdist)
    # Go from black-and-white to defined colours. Flip axes to image-like.
    map(res) do pix
        pix == true && return RGBA{N0f8}(0.714, 0.333, 0.0, 1)
        RGBA{N0f8}(0., 0, 0, 0)
    end
end
function _elev_contours!(res::T1, zs::T2, 
        bbuf::T1, minlen::Int64, vthick::T3, elevation_spacings::T3) where {T1 <: Matrix{Gray{Bool}},
            T2 <: Matrix{RGB{Float32}},
            T3 <: Vector{Int64}}
    for (t, Δz) in zip(vthick, elevation_spacings)
        # Overwrite bbuf with pixels on contours.
        mapwindow!(func_elev_contour(Δz), bbuf, zs, (3, 3))
        # We have tried to reduce contours on bumps. Now remove most of the rest:
        # Our very clever thinning, despeckling and thickening
        # (and this part takes ~85% of the time in this loop):
        bbuf .= strokeify(bbuf, t, minlen)
        # Treat 'false' as transparent. Overlay bbuf on res and write in-place to res.
        map!(BlendLighten, res, res, bbuf)
    end
    res
end

function terrain_gradient(zs::Matrix{Float32})
    # Prepare a wide low-pass filter. Window size:
    w = 49
    # Coefficients, including a window.
    c = Float32.(fir_lp_coefficients(w))
    # Elevations smoothed by lp-filter
    zf = imfilter(zs, (c, transpose(c)), FIRTiled());
    # Find nearby directional steepness, i.e. not the steepness of the
    # local stone, tree or shed.
    imgradients(zf, KernelFactors.ando5, "replicate")
end

"""
    strokeify(bw, thickness, minlen)

Applies a series of filters that change a collection of adjacent pixels to
'pen strokes' with the given thickness in pixels. 

Increase minlen to drop 'dashed elevation contours' in forests, built-up areas.
"""
function strokeify(bw, thickness, minlen)
    thickness == 0 && return bw
    if (thickness - 1 ) % 2 !== 0
        throw(ArgumentError("thickness = $thickness ∉ [0, 1, 3, ...] "))
    end
    # Despeckle
    remove_isolated_pixels!(bw)
    # Reduce all the 'true' regions to one pixel width.
    stroked = reinterpret(Gray{Bool}, thinning(reinterpret(Bool, bw), algo = GuoAlgo()))
    # Remove too short contours, (typically, minlen = 5)
    remove_small_islands!(stroked, minlen)
    # Make all the thin regions as thick as specified
    thickness == 1 && return stroked
    dilate(stroked; r = thickness ÷ 2)
end

"""
    func_elev_contour(Δc::Float32)
    func_elev_contour(Δc::Float32)
    ---> funtion

Δc is the vertical distance between contours, typically 20, 100 or 1000 m.
"""
func_elev_contour(Δc) = func_elev_contour(Float32(Δc))
function func_elev_contour(Δc::Float32)
    max_from_contour = 3.0f0
    isnot =  Gray{Bool}(false)
    is = Gray{Bool}(true)
    # 
    f = let Δc = Δc, max_from_contour = max_from_contour, isnot =  isnot, is = is
        (M) -> let
            @assert size(M) == (3, 3)
            # If rows in M correspond to south -> north
            # and cols in M correspond to west -> east
            # _ n _
            # w z e  
            # _ s _
            _, w, _, n, z, s, _, e, _ = red.(M)
            # No elevation contour line at sea level
            z < Δc && return isnot
            # Elevation at centre relative to the nearest contour
            # Δz will vary as a sawtooth function between (-Δc / 2), below the contour,
            # and Δc / 2  above the nearest contour.
            Δz = mod(z - Δc / 2, Δc) - Δc / 2
            # Are we even roughly close to the contour? If not, we can decide fast.
            abs(Δz) > max_from_contour  && return isnot
            # Are any of the neighbours not on the same 'sawtooth'? 
            # That might lead to artifacts, and we better judge that this is no contour.
            # This will reduce the overlapping contour pixels in very steep parts.
            abs(z - e) >= Δc / 2 && return isnot
            abs(z - w) >= Δc / 2 && return isnot
            abs(z - n) >= Δc / 2 && return isnot
            abs(z - s) >= Δc / 2 && return isnot
              # Elevation at immediate neighbours relative to the nearest contour
            Δw = mod(w - Δc / 2, Δc) - Δc / 2
            Δn = mod(n - Δc / 2, Δc) - Δc / 2
            Δs = mod(s - Δc / 2, Δc) - Δc / 2
            Δe = mod(e - Δc / 2, Δc) - Δc / 2
            # If west-east elevations do not cross the contour,
            # then sign(Δw) * sign(Δe) == 1
            cross_we = sign(Δw) * sign(Δe)
            # If south-north elevations do not cross the contour,
            # then sign(Δs) * sign(Δn) == 1
            cross_sn = sign(Δs) * sign(Δn)
            # If both are non-crossing, this is no contour, for sure
            cross_we == cross_sn == 1 && return isnot
            # We still would have too many false positives, i.e. 
            # 'noise', mostly due to trees and houses.
            # So also consider with regional gradients. Values 
            # are taken from highly smoothed gradient, at the 
            # centre pixel):
            pix = M[5]
            pix_y´ = green(pix)
            pix_x´ = blue(pix)
            # Remember 'y' increases downwards,
            # so: sign(wy´) == 1 slopes upward in a southerly direction. 
            # Let's eliminate local upcrossings in regions of downslopes
            sign(pix_y´) !== sign(s - n) && return isnot
            sign(pix_x´) !== sign(e - w) && return isnot
            # This is, as far as we can tell, a valid contour crossing.
            return is
        end
    end
    f
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
