using Test
using BitmapMaps
using BitmapMaps: SheetMatrixBuilder

# SheetMatrixBuilder corresponds to resource/matrix_sheet_cell_utm.svg
bm_cell_width, bm_cell_height, sheet_cell_width, sheet_cell_height, sw_corner, cell_to_utm_factor = (9, 8, 3, 4, (44000, 6909047), 2)
smb = SheetMatrixBuilder(bm_cell_width, bm_cell_height, sheet_cell_width, sheet_cell_height, sw_corner, cell_to_utm_factor, 191, "BitmapMaps/test")
@test length(smb) == 6


@test northwest_corner(smb) .- sw_corner ==  (0, bm_cell_height) .* cell_to_utm_factor
@test southwest_corner(smb) .- sw_corner ==  (0, 0) 
@test northeast_corner(smb) .- sw_corner ==  (bm_cell_width, bm_cell_height) .* cell_to_utm_factor
@test southeast_corner(smb) .- sw_corner ==  (bm_cell_width, 0) .* cell_to_utm_factor


p = smb[2, 1]
@test northwest_corner(p) .- sw_corner == (0, 2 * sheet_cell_height) .* cell_to_utm_factor
@test southwest_external_corner(p) .- sw_corner == (0, sheet_cell_height) .* cell_to_utm_factor
@test northeast_external_corner(p) .- sw_corner == (sheet_cell_width, 2 * sheet_cell_height) .* cell_to_utm_factor
@test southeast_external_corner(p) .- sw_corner == (sheet_cell_width, sheet_cell_height) .* cell_to_utm_factor

@test southwest_internal_corner(p) .- sw_corner == (0, sheet_cell_height + 1) .* cell_to_utm_factor
@test northeast_internal_corner(p) .- sw_corner == (sheet_cell_width - 1, sheet_cell_height * 2) .* cell_to_utm_factor
@test southeast_internal_corner(p) .- sw_corner == (sheet_cell_width - 1, sheet_cell_height + 1) .* cell_to_utm_factor

@test geo_centre(smb) isa Tuple{Float64, Float64}
@test geo_centre(smb[1,1]) isa Tuple{Float64, Float64}
@test geo_grid_centre_single(smb) isa Tuple{Int64, Int64}
@test geo_grid_centre_single(smb[1,1]) isa Tuple{Int64, Int64}

@test geo_centre(smb) == (44009, 6909055)
@test geo_grid_centre_single(smb) == (44009, 6909056) 
@test geo_centre(smb[2, 2]) == (44009, 6909059)
@test geo_grid_centre_single(smb[2,2]) == (44009, 6909060) 

@test all(geo_centre(smb)                  .> southwest_corner(smb))
@test all(geo_centre(smb[1,1])             .> southwest_internal_corner(smb[1,1]))
@test all(geo_grid_centre_single(smb)      .> southwest_corner(smb))
@test all(geo_grid_centre_single(smb[1,1]) .> southwest_internal_corner(smb[1,1]))

@test all(geo_centre(smb)                  .< northeast_internal_corner(smb))
@test all(geo_centre(smb[1,1])             .< northeast_internal_corner(smb[1,1]))
@test all(geo_grid_centre_single(smb)      .< northeast_internal_corner(smb))
@test all(geo_grid_centre_single(smb[1,1]) .< northeast_internal_corner(smb[1,1]))


# Check that a sheet's folder name matches its northeast_external_corner
@test tryparse.(Int, split(smb[2,2].pthsh, ' ')[end - 1: end]) == collect(northeast_external_corner(smb[2, 2]))


@test geo_area(smb) == smb.bm_cell_width * smb. bm_cell_height * cell_to_utm_factor^2
@test geo_area(smb[1]) == (smb.bm_cell_width * smb. bm_cell_height * cell_to_utm_factor^2) / (smb.nrows * smb.ncols)

# TODO bad name. We don't know the values, nonzero doesn't make sense. Make 'raster-string' or 'bb-string'?
@test nonzero_raster_string(smb) == "(44000 6909047)-(44018 6909063)"
@test nonzero_raster_string(smb[1,1]) == "(44000 6909047)-(44006 6909055)"
@test replace(polygon_string(smb), "\t" => "") == """
    POLYGON ((44000 6909047, 44006 6909047, 44006 6909055, 44000 6909055, 44000 6909047),
    (44000 6909055, 44006 6909055, 44006 6909063, 44000 6909063, 44000 6909055),
    (44006 6909047, 44012 6909047, 44012 6909055, 44006 6909055, 44006 6909047),
    (44006 6909055, 44012 6909055, 44012 6909063, 44006 6909063, 44006 6909055),
    (44012 6909047, 44018 6909047, 44018 6909055, 44012 6909055, 44012 6909047),
    (44012 6909055, 44018 6909055, 44018 6909063, 44012 6909063, 44012 6909055))"""

show_augmented_properties(smb)
show_augmented_properties(smb[end, end])