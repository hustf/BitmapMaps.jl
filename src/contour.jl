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
    # Save
    ffna = joinpath(fofo, CONTOUR_FNAM)
    @debug "    Saving $ffna"
    save_png_with_phys(ffna, res)
    true
end
function _elev_contours(fofo, cell_iter, cell2utm, minlen, vthick, vdist)
    # Read elevation
    g = readclose(joinpath(fofo, CONSOLIDATED_FNAM))
    #
    # Prepare iterators and buffers, reduce source resolution
    #
    ny, nx = size(cell_iter)
    source_indices = (1:cell2utm:(ny  * cell2utm), 1:cell2utm:(nx * cell2utm))
    si = CartesianIndices(source_indices)
    # We store three arrays in an 'image' in order to use the fast image filtering functions later.
    #    red:   is_forest
    #    green: unfiltered elevation
    #    blue:  highly smoothed elevation
    # Note that although we use the full resolution for identifying forest,
    # we reduce the resolution in the three-channel image.
    img = RGB{Float32}.(Float32.(is_forest(g)[si]), # red
            transpose(g.A[:, :, 1])[si],           # green
            Float32.(imfilter(transpose(g.A[:, :, 1])[si], Kernel.gaussian(49)))) # blue
    #rm = maximum(red.(img))
    #gm = maximum(green.(img))
    #bm = maximum(blue.(img))
    #img1 = RGB{N0f8}.(red.(img) ./ rm, green.(img) ./ gm, blue.(img) ./ bm)
    #display(img1)
    #ffna = joinpath(fofo, "temp.png")
    #@debug "    Saving $ffna"
    #save_png_with_phys(ffna, img1)
    res = __elev_contours(img, minlen, vthick, vdist)
    # Go from black-and-white to defined colours. Flip axes to image-like.
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
    # Despeckle
    # TEMP remove_isolated_pixels!(bw)
    # Reduce regions to one pixel width.
    stroked = reinterpret(Gray{Bool}, thinning(reinterpret(Bool, bw), algo = GuoAlgo()))
    # Remove too short contours, (typically, minlen = 5)
    # TEMP remove_small_islands!(stroked, minlen)
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
            @assert size(M) == (3, 3)
            # If rows in M correspond to south -> north
            # and cols in M correspond to west -> east
            # _ n _
            # w z e  
            # _ s _
            is_forest = red(M[5]) == 1.0f0 
            if is_forest 
                # Use smoothed terrain values.
                # We also compensate for a constant, assumed mean forest height.
                # This introduces a error at the edge of forests,
                # which could be masked by fading in the subtraction (via segmented distance calculations).
                nw, w, sw, n, z, s, ne, e, se = blue.(M) .- 4f0
            else
                nw, w, sw, n, z, s, ne, e, se = green.(M)      # Exact terrain
            end
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
