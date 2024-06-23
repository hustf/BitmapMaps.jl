# Step in pipeline.
# Creates utm grid lines overlay.
# # Output is an image file per sheet, for manual touch-up.

"""
    grid_overlay(sb::SheetBuilder)
    grid_overlay(fofo, cell_iter, cell_to_utm_factor)

"""
grid_overlay(sb::SheetBuilder) = grid_overlay(full_folder_path(sb), sb.cell_iter, cell_to_utm_factor(sb), sb.f_I_to_utm)
function grid_overlay(fofo, cell_iter, cell2utm, f_I_to_utm)
    if isfile(joinpath(fofo, COMPOSITE_FNAM))
        @debug "$COMPOSITE_FNAM in $fofo already exists. Exiting `grid_overlay`."
        return true
    end
    if isfile(joinpath(fofo, GRID_FNAM))
        @debug "$GRID_FNAM in $fofo already exists. Exiting `grid_overlay`."
        return true
    end
    res = _grid_utm(fofo, cell_iter, cell2utm, f_I_to_utm)
    ffna = joinpath(fofo, GRID_FNAM)
    # We won't ever print this. The value won't be used. So we specify a standard 300 dpi, disregarding user specs
    # for the bitmapmap
    density_pt_m⁻¹ = 11811
    @debug "Saving $ffna"
    save_png_with_phys(ffna, res, density_pt_m⁻¹)
    true
end
function _grid_utm(fofo, cell_iter, cell2utm, f_I_to_utm)
    # Establish output image matrix
    ny, nx = size(cell_iter)
    linecol = RGBA{N0f8}(1.0, 0.843, 0.0, 0.651)
    transpcol = RGBA{N0f8}(0, 0, 0, 0)
    @debug "Render grid lines"
    fg = generate_grid_utm_func(linecol, transpcol, f_I_to_utm)
    grid = map(fg, cell_iter)
    display(grid)
    grid
end

function generate_grid_utm_func(linecol, transpcol, f_I_to_utm; width_1000 = 8)
    (I::CartesianIndex) -> begin
        easting, northing = f_I_to_utm(I)
        # Wide contour line at 1000m, 2000m
        easting % 1000 < (width_1000 ) && return linecol
        northing % 1000 < (width_1000 ) && return linecol
        transpcol
    end
end
