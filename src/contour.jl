# Step in pipeline.
# Creates elevation contour lines overlay from elevation data.
# Output is an image file per sheet, for manual touch-up.
# Applies smoothing to terrain where that is needed, but
# keeps details where there is no forest, buildings and the like.

"""
    contour_lines_overlay(sb::SheetBuilder)
    contour_lines_overlay(fofo, cell_iter, cell_to_utm_factor)
    ---> Bool

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
    # Save
    ffna = joinpath(fofo, CONTOUR_FNAM)
    @debug "    Saving $ffna"
    save_png_with_phys(ffna, res)
    true
end
function _elev_contours(fofo, cell_iter, cell2utm, minlen, vthick, vdist)
    # Read elevation
    z = elevation_full_res(fofo)
    #
    # Prepare iterators and buffers, reduce source resolution
    #
    ny, nx = size(cell_iter)
    source_indices = (1:cell2utm:(ny  * cell2utm), 1:cell2utm:(nx * cell2utm))
    si = CartesianIndices(source_indices)
    # We store three arrays in an 'image' in order to use the fast image filtering functions later.
    #    red:   bumpy_patch
    rr = Float32.(bumpy_patch(z, si))
    #    green: unfiltered elevation
    gg = z[si]
    #    blue:  highly smoothed elevation
    bb = Float32.(imfilter(z, Kernel.gaussian(49))[si])
    # Note that although we use the full resolution for identifying bumpy patches,
    # we reduce the resolution in this three-channel image:
    img = RGB{Float32}.(rr, gg,  bb) # blue
    # Boolean countours image:
    res = __elev_contours(img, minlen, vthick, vdist)
    # An alternative topo colour which would show better on green:
    #     RGB 0.71, 0.572, 0. 157
    # We consider topo lines are unimportant in such places, so better to use one colour.
    #
    # Go from black-and-white to defined colours.
    map(res) do pix
        pix == true && return RGBA{N0f8}(0.714, 0.333, 0.0, 1)
        RGBA{N0f8}(0., 0, 0, 0)
    end
end
function __elev_contours(img, minlen, vthick, vdist)
    # Pre-allocate boolean buffer (we make one countour distance at a time)
    bbuf = Array{Gray{Bool}}(undef, size(img)...)
    # Pre-allocate boolean result image
    res = zeros(Gray{Bool}, size(img)...)
    # With all buffers ready, call the inner function
    _elev_contours!(res, img, bbuf, minlen, vthick, vdist)
    res
end
function _elev_contours!(res::T1, img::T2,
    bbuf::T1, minlen::Int64, vthick::T3, elevation_spacings::T3) where {T1 <: Matrix{Gray{Bool}},
        T2 <: Matrix{RGB{Float32}},
        T3 <: Vector{Int64}}
    for (t, Δz) in zip(vthick, elevation_spacings)
        # Overwrite bbuf with pixels on contours.
        mapwindow!(func_elev_contour(Δz), bbuf, img, (3, 3))
        # We have tried to reduce contours on bumps. Now remove most of the rest:
        # Our very clever thinning, despeckling and thickening
        # (and this part takes ~85% of the time in this loop):
        bbuf .= strokeify(bbuf, t, minlen)
        # Treat 'false' as transparent. Overlay bbuf on res and write in-place to res.
        map!(BlendLighten, res, res, bbuf)
    end
    res
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
    # Reduce regions to one pixel width.
    stroked = reinterpret(Gray{Bool}, thinning(reinterpret(Bool, bw), algo = GuoAlgo()))
    # Remove too short contours, (typically, minlen = 5)
    remove_small_strokes!(stroked, minlen)
    # Make all the thin regions as thick as specified
    thickness == 1 && return stroked
    dilate(stroked; r = thickness ÷ 2)
end

"""
    func_elev_contour(Δc::Float32)
    ---> function

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
            #    red:   bumpy_patch
            #    green: unfiltered elevation
            #    blue:  highly smoothed elevation
            @assert size(M) == (3, 3)
            # If rows in M correspond to south -> north
            # and cols in M correspond to west -> east
            # _ n _
            # w z e
            # _ s _
            #
            # If all pixels are in the forest, we apply the blue channel fully.
            # Otherwise, we interpolate according to 'forest factor':
            forest_factor = sum(x -> x > 0, red.(M)) / 9
            # Interpolate from forest_factor.
            # Note that we subtract 4 m from the elevation in bumpy patches. This is
            # often too much in upper elevations, and too little in lower elevations
            # with tall forest and houses. We aren't able to estimate the actual
            # height of the forest.
            nw, w, sw, n, z, s, ne, e, se  = forest_factor .* (blue.(M) .- 4f0 ) .+ (1 - forest_factor)  .* green.(M)
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
            # Elevation at immediate neighbours relative to the nearest contour.
            # If e.g. west-east elevations cross the contour,
            # then sign(Δw) * sign(Δe) == - 1
            Δw  = mod(w - Δc / 2, Δc) - Δc / 2
            Δe  = mod(e - Δc / 2, Δc) - Δc / 2
            sign(Δw) * sign(Δe) == -1 && return is
            Δn  = mod(n - Δc / 2, Δc) - Δc / 2
            Δs  = mod(s - Δc / 2, Δc) - Δc / 2
            sign(Δs) * sign(Δn) == -1 && return is
            Δnw = mod(nw - Δc / 2, Δc) - Δc / 2
            Δse = mod(se - Δc / 2, Δc) - Δc / 2
            sign(Δnw) * sign(Δse) == -1 && return is
            Δsw = mod(sw - Δc / 2, Δc) - Δc / 2
            Δne = mod(ne - Δc / 2, Δc) - Δc / 2
            sign(Δsw) * sign(Δne) == -1 && return is
            return isnot
        end
    end
    f
end


function remove_small_strokes!(img, max_pixels)
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