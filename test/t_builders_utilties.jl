using Test
using BitmapMaps
using BitmapMaps: SheetMatrixBuilder

# Corresponds to resource/matrix_sheet_cell_utm.svg
sw_corner, nrc, cell_to_utm_factor = (44000, 6909047), (2,3), 2
sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, pth = 3, 4, 1000, "nopath"
smb = SheetMatrixBuilder(sw_corner, nrc, cell_to_utm_factor, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, pth)

@test BitmapMaps.nrows(smb) == 2
@test BitmapMaps.ncols(smb) == 3


@test northwest_corner(smb) .- sw_corner ==  (0, BitmapMaps.bm_height_cell(smb)) .* cell_to_utm_factor
@test southwest_corner(smb) .- sw_corner ==  (0, 0) 
@test northeast_corner(smb) .- sw_corner ==  (BitmapMaps.bm_width_cell(smb), BitmapMaps.bm_height_cell(smb)) .* cell_to_utm_factor
@test southeast_corner(smb) .- sw_corner ==  (BitmapMaps.bm_width_cell(smb), 0) .* cell_to_utm_factor


p = smb[2, 1]
@test northwest_corner(p) .- sw_corner == (0, 2 * BitmapMaps.sheet_height_cell(smb)) .* cell_to_utm_factor
@test southwest_external_corner(p) .- sw_corner == (0, BitmapMaps.sheet_height_cell(smb)) .* cell_to_utm_factor
@test northeast_external_corner(p) .- sw_corner == (BitmapMaps.sheet_width_cell(smb), 2 * BitmapMaps.sheet_height_cell(smb)) .* cell_to_utm_factor
@test southeast_external_corner(p) .- sw_corner == (BitmapMaps.sheet_width_cell(smb), BitmapMaps.sheet_height_cell(smb)) .* cell_to_utm_factor

@test southwest_internal_corner(p) .- sw_corner == (0, BitmapMaps.sheet_height_cell(smb) + 1) .* cell_to_utm_factor
@test northeast_internal_corner(p) .- sw_corner == (BitmapMaps.sheet_width_cell(smb) - 1, BitmapMaps.sheet_height_cell(smb) * 2) .* cell_to_utm_factor
@test southeast_internal_corner(p) .- sw_corner == (BitmapMaps.sheet_width_cell(smb) - 1, BitmapMaps.sheet_height_cell(smb) + 1) .* cell_to_utm_factor

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


@test geo_area(smb) == BitmapMaps.bm_width_cell(smb) * BitmapMaps.bm_height_cell(smb) * cell_to_utm_factor^2
@test geo_area(smb[1]) == (BitmapMaps.bm_width_cell(smb) * BitmapMaps.bm_height_cell(smb) * cell_to_utm_factor^2) / (BitmapMaps.nrows(smb) * BitmapMaps.ncols(smb))

@test bbox_external_string(smb) == "(44000 6909047)-(44018 6909063)"
@test bbox_external_string(smb[1,1]) == "(44000 6909047)-(44006 6909055)"
@test replace(polygon_string(smb), "\t" => "") == """
    POLYGON ((44000 6909047, 44006 6909047, 44006 6909055, 44000 6909055, 44000 6909047),
    (44000 6909055, 44006 6909055, 44006 6909063, 44000 6909063, 44000 6909055),
    (44006 6909047, 44012 6909047, 44012 6909055, 44006 6909055, 44006 6909047),
    (44006 6909055, 44012 6909055, 44012 6909063, 44006 6909063, 44006 6909055),
    (44012 6909047, 44018 6909047, 44018 6909055, 44012 6909055, 44012 6909047),
    (44012 6909055, 44018 6909055, 44018 6909063, 44012 6909063, 44012 6909055))"""

show_derived_properties(smb)
show_derived_properties(smb[end, end])