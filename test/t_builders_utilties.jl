using Test
using BitmapMaps
using BitmapMaps: SheetMatrixBuilder

# SheetMatrixBuilder corresponds to resource/matrix_sheet_pix_utm.svg
bm_pix_width, bm_pix_height, sheet_pix_width, sheet_pix_height, sw_corner, pix_to_utm_factor = (9, 8, 3, 4, (43999, 6909048), 2)
smb = SheetMatrixBuilder(bm_pix_width, bm_pix_height, sheet_pix_width, sheet_pix_height, sw_corner, pix_to_utm_factor, "test")
@test length(smb) == 6


@test northwest_corner(smb) .- sw_corner ==  (0, bm_pix_height) .* pix_to_utm_factor
@test southwest_corner(smb) .- sw_corner ==  (0, 0) 
@test northeast_corner(smb) .- sw_corner ==  (bm_pix_width, bm_pix_height) .* pix_to_utm_factor
@test southeast_corner(smb) .- sw_corner ==  (bm_pix_width, 0) .* pix_to_utm_factor


p = smb[2, 1]
@test northwest_corner(p) .- sw_corner == (0, 2 * sheet_pix_height) .* pix_to_utm_factor
@test southwest_external_corner(p) .- sw_corner == (0, sheet_pix_height) .* pix_to_utm_factor
@test northeast_external_corner(p) .- sw_corner == (sheet_pix_width, 2 * sheet_pix_height) .* pix_to_utm_factor
@test southeast_external_corner(p) .- sw_corner == (sheet_pix_width, sheet_pix_height) .* pix_to_utm_factor

@test southwest_internal_corner(p) .- sw_corner == (0, sheet_pix_height + 1) .* pix_to_utm_factor
@test northeast_internal_corner(p) .- sw_corner == (sheet_pix_width - 1, sheet_pix_height * 2) .* pix_to_utm_factor
@test southeast_internal_corner(p) .- sw_corner == (sheet_pix_width - 1, sheet_pix_height + 1) .* pix_to_utm_factor

@test geo_centre(smb) isa Tuple{Float64, Float64}
@test geo_centre(smb[1,1]) isa Tuple{Float64, Float64}
@test geo_grid_centre_single(smb) isa Tuple{Int64, Int64}
@test geo_grid_centre_single(smb[1,1]) isa Tuple{Int64, Int64}

@test geo_centre(smb) == (44008, 6909056)
@test geo_grid_centre_single(smb) == (44008, 6909057) 
@test geo_centre(smb[2, 2]) == (44008, 6909060)
@test geo_grid_centre_single(smb[2,2]) == (44008, 6909061) 

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


@test geo_area(smb) == smb.bm_pix_width * smb. bm_pix_height * pix_to_utm_factor^2
@test geo_area(smb[1]) == (smb.bm_pix_width * smb. bm_pix_height * pix_to_utm_factor^2) / (smb.nrows * smb.ncols)

@test bounding_box_external_string(smb) == "(43999 6909048)-(44017 6909064)"
@test bounding_box_external_string(smb[1,1]) == "(43999 6909048)-(44005 6909056)"
@test replace(BitmapMaps.bounding_box_polygon_string(smb), "\t" => "") == """
    POLYGON ((43999 6909048, 44005 6909048, 44005 6909056, 43999 6909056, 43999 6909048),
    (43999 6909056, 44005 6909056, 44005 6909064, 43999 6909064, 43999 6909056),
    (44005 6909048, 44011 6909048, 44011 6909056, 44005 6909056, 44005 6909048),
    (44005 6909056, 44011 6909056, 44011 6909064, 44005 6909064, 44005 6909056),
    (44011 6909048, 44017 6909048, 44017 6909056, 44011 6909056, 44011 6909048),
    (44011 6909056, 44017 6909056, 44017 6909064, 44011 6909064, 44011 6909056))"""

show_augmented_properties(smb)
show_augmented_properties(smb[end, end])