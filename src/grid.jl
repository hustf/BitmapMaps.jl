# Step in pipeline.
# Creates utm grid lines overlay.
# Output is an image file per sheet, for manual touch-up.

"""
    grid_overlay(sb::SheetBuilder)
    grid_overlay(fofo, cell_iter, cell_to_utm_factor)

"""
function grid_overlay(sb::SheetBuilder)
    thick_grid = cell_to_utm_factor(sb) * get_config_value("UTM grid", "Grid line thickness", Int)
    spacing_grid = get_config_value("UTM grid", "Grid spacing [m]", Int)
    zone_no = get_config_value("UTM grid", "Zone for elevation data", Int)
    # And maybe one for wgs84?? 
    grid_overlay(full_folder_path(sb), sb.cell_iter, sb.f_I_to_utm, thick_grid, spacing_grid, zone_no)
end
function grid_overlay(fofo, cell_iter, f_I_to_utm, thick_grid, spacing_grid, zone_no)
    if isfile(joinpath(fofo, GRID_FNAM))
        @debug "    $GRID_FNAM in $fofo \n           already exists. Exiting `grid_overlay`"
        return true
    end
    res = _grid_utm(cell_iter, f_I_to_utm, thick_grid, spacing_grid, zone_no)
    # Save
    ffna = joinpath(fofo, GRID_FNAM)
    @debug "    Saving $ffna"
    save_png_with_phys(ffna, res)
    true
end

function _grid_utm(cell_iter, f_I_to_utm, thick_grid, spacing_grid, zone_no)
    linecol = RGBA{N0f8}(1.0, 0.843, 0.0, 0.651)
    transpcol = RGBA{N0f8}(0, 0, 0, 0)
    @debug "    Render grid lines"
    fg = func_grid_utm(linecol, transpcol, f_I_to_utm, thick_grid, spacing_grid, zone_no)
    # Slow, but no matter
    map(fg, cell_iter)
end

function func_grid_utm(linecol, transpcol, f_I_to_utm, thick_grid, spacing_grid, zone_no)
    # This function is a CoordinateTransformations.ComposedTransformation.
    # It could be inverted for faster (more complicated) code.
    # This ineffective approach takes ~9 seconds per A4 sheet, 
    # but seldom needs repetition. This is neglectible compared to 
    # finding water surfaces.
    transform = UTMZfromLLA(wgs84) ∘ LLAfromUTMZ(wgs84)
    function to_localgrid(I)
        a = UTMZ(UTM(f_I_to_utm(I)...), zone_no, true)
        b = transform(a)
        Int(round(b.x)), Int(round(b.y))
    end
    lower_limit = div(thick_grid, 2) - 1
    upper_limit = spacing_grid - div(thick_grid, 2)
    (I::CartesianIndex) -> begin
        easting, northing = to_localgrid(I)
        remainder = easting % spacing_grid
        remainder <= lower_limit || remainder >= upper_limit && return linecol
        remainder = northing % spacing_grid
        remainder <= lower_limit || remainder >= upper_limit && return linecol
        transpcol
    end
end
