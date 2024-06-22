# Creates elevation contour lines overlay from elevation data.
# # Output is an image file per sheet, for manual touch-up.

"""
    contour_lines_overlay(sb::SheetBuilder)
    contour_lines_overlay(fofo, cell_iter, cell_to_utm_factor)

"""
contour_lines_overlay(sb::SheetBuilder) = contour_lines_overlay(full_folder_path(sb), sb.cell_iter, cell_to_utm_factor(sb))
function contour_lines_overlay(fofo, cell_iter, cell2utm)
    if isfile(joinpath(fofo, CONTOUR_FNAM))
        @debug "$CONTOUR_FNAM in $fofo already exists. Exiting `contour_lines_overlay`."
        return true
    end
    _elev_contour(fofo, cell_iter, cell2utm)
end
function _elev_contour(fofo, cell_iter, cell2utm)
    # Read elevation
    g = readclose(joinpath(fofo, CONSOLIDATED_FNAM))
    za = g.A[:, :, 1]
    # Establish output image matrix
    ny, nx = size(cell_iter)
    source_indices = (1:cell2utm:(nx * cell2utm), 1:cell2utm:(ny * cell2utm))
    linecol = RGBA{N0f8}(0.173, 0.192, 0.255, 1)
    transpcol = RGBA{N0f8}(0, 0, 0, 0)
    fc = generate_elev_contour_func(linecol, transpcol)
    @debug "Render elevation contour lines"
    contour = mapwindow(fc, za, (9, 9), indices = source_indices)
    display(transpose(contour))
    ffna = joinpath(fofo, CONTOUR_FNAM)
    # We won't ever print this. The value won't be used. So we specify a standard 300 dpi, disregarding user specs
    # for the bitmapmap
    density_pt_m⁻¹ = 11811
    @debug "Saving $ffna"
    save_png_with_phys(ffna, transpose(contour), density_pt_m⁻¹)
    true
end

function generate_elev_contour_func(linecol, transpcol; width_1000 = 8, width_100 = 4)
    (M::Matrix) -> @inbounds begin
        # This function takes elevation in a local grid with 1m spacing, regardless
        # of the output image density
        @assert size(M) == (9, 9)
        # We consider the slope because we want equally
        # wide contour lines.
        # The slope varies somewhat too much over a distance of 2 m.
        # Hence, we use a wide span (a large window)
        span = 8
        w = M[1, 5]
        e = M[9, 5]
        n = M[5, 1]
        z = M[5, 5]
        s = M[5, 9]
        deriv_south_north = (n - s) / span
        deriv_east_west = (e - w) / span
        gradient_magn = hypot(deriv_south_north, deriv_east_west)
        # No elevation contour line at sea level 
        z < 90 && return transpcol
        # Wide contour line at 1000m, 2000m
        z % 1000 < (width_1000 * gradient_magn) && return linecol
        # Narrow contour line at 100m, 200m
        z % 100 < (width_100 * gradient_magn) && return linecol
        transpcol
    end
end
