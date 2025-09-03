using Test
using BitmapMaps
using BitmapMaps: cell_to_utm_factor, width_cell, height_cell
using BitmapMaps: func_utm_to_sheet_index, func_utm_to_cell_index, geo_width_height

# Corresponds to resource/matrix_sheet_cell_utm.svg
sw_corner, nrc, cell2utm = (44000, 6909047), (2,3), 2
sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, pth = 3, 4, 1000, "nopath"
smb = SheetMatrixBuilder(sw_corner, nrc, cell2utm, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, pth)

@test BitmapMaps.nrows(smb) == 2
@test BitmapMaps.ncols(smb) == 3


@test northwest_corner(smb) .- sw_corner ==  (0, BitmapMaps.height_cell(smb)) .* cell2utm
@test southwest_corner(smb) .- sw_corner ==  (0, 0)
@test northeast_corner(smb) .- sw_corner ==  (BitmapMaps.width_cell(smb), BitmapMaps.height_cell(smb)) .* cell2utm
@test southeast_corner(smb) .- sw_corner ==  (BitmapMaps.width_cell(smb), 0) .* cell2utm


p = smb[2, 1]
@test northwest_corner(p) .- sw_corner == (0, 2 * BitmapMaps.sheet_height_cell(smb)) .* cell2utm
@test southwest_external_corner(p) .- sw_corner == (0, BitmapMaps.sheet_height_cell(smb)) .* cell2utm
@test northeast_external_corner(p) .- sw_corner == (BitmapMaps.sheet_width_cell(smb), 2 * BitmapMaps.sheet_height_cell(smb)) .* cell2utm
@test southeast_external_corner(p) .- sw_corner == (BitmapMaps.sheet_width_cell(smb), BitmapMaps.sheet_height_cell(smb)) .* cell2utm
@test southwest_internal_corner(p) .- sw_corner == (0, BitmapMaps.sheet_height_cell(smb) + 1) .* cell2utm
@test northeast_internal_corner(p) .- sw_corner == (BitmapMaps.sheet_width_cell(smb) - 1, BitmapMaps.sheet_height_cell(smb) * 2) .* cell2utm
@test southeast_internal_corner(p) .- sw_corner == (BitmapMaps.sheet_width_cell(smb) - 1, BitmapMaps.sheet_height_cell(smb) + 1) .* cell2utm

p = smb[1, 2]
wp = width_cell(p) * cell_to_utm_factor(p)
hp = height_cell(p) * cell_to_utm_factor(p)
Δ = cell_to_utm_factor(p)
@test northwest_corner(p) .- sw_corner == (wp, hp)
@test southwest_external_corner(p) .- sw_corner == ( wp, 0)
@test northeast_external_corner(p) .- sw_corner == (2wp, hp)
@test southeast_external_corner(p) .- sw_corner == (2wp, 0)
@test southwest_internal_corner(p) .- sw_corner == ( wp, Δ)
@test northeast_internal_corner(p) .- sw_corner == (2wp - Δ, hp)
@test southeast_internal_corner(p) .- sw_corner == (2wp - Δ,  Δ)




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
@test tryparse.(Int, split(split(smb[2,2].pthsh, "__")[end], '-')) == collect(northeast_external_corner(smb[2, 2]))


@test geo_area(smb) == BitmapMaps.width_cell(smb) * BitmapMaps.height_cell(smb) * cell2utm^2
@test geo_area(smb[1]) == (BitmapMaps.width_cell(smb) * BitmapMaps.height_cell(smb) * cell2utm^2) / (BitmapMaps.nrows(smb) * BitmapMaps.ncols(smb))

@test bbox_external_string(smb) == "(44000 6909049)-(44016 6909063)"
@test bbox_external_string(smb[1,1]) == "(44000 6909049)-(44004 6909055)"
@test polygon_string(smb[1,2]) == "MULTIPOLYGON (\n                   ((44006 6909047, 44012 6909047, 44012 6909055, 44006 6909055, 44006 6909047)))"
@test polygon_string(smb) == "MULTIPOLYGON (\n                   ((44000 6909047, 44006 6909047, 44006 6909055, 44000 6909055, 44000 6909047)),\n                   ((44000 6909055, 44006 6909055, 44006 6909063, 44000 6909063, 44000 6909055)),\n                   ((44006 6909047, 44012 6909047, 44012 6909055, 44006 6909055, 44006 6909047)),\n                   ((44006 6909055, 44012 6909055, 44012 6909063, 44006 6909063, 44006 6909055)),\n                   ((44012 6909047, 44018 6909047, 44018 6909055, 44012 6909055, 44012 6909047)),\n                   ((44012 6909055, 44018 6909055, 44018 6909063, 44012 6909063, 44012 6909055)))"


##############################################################
# Test conversion from geographical coordinates to sheet index
##############################################################

w_utm, h_utm = geo_width_height(smb)
ws_utm, hs_utm = geo_width_height(smb[1])
nw_utm = northwest_corner(smb[1])
# The coordinates defined in the figure resource/matrix_sheet_cell_utm
e, n = southwest_corner(smb)


# Generated utm to sheet index function
f = func_utm_to_sheet_index(smb)
fs(e, n) = f((e, n))
# The southwest corner, touched by four sheets, belongs to a sheet on the south side of the sheet matrix.
# Whose coordinates are outside the sheet matrix.
@test_throws ArgumentError fs(e, n) 
ei, ni = southwest_internal_corner(smb)
@test fs(ei, ni) == (1, 1)
# The sheet row is 1 for this range of utm northing.
@test [fs(e, n + i)[1] for i in 1:hs_utm] == ones(Int, hs_utm)
# The sheet row is 2 for this range of utm northing.
@test [fs(e, n + i)[1] for i in (1+hs_utm):2hs_utm] == 2 .* ones(Int, hs_utm)
# The northheast corner, touched by four sheets, belongs to the sheet on the east side of the sheet matrix.
# Whose coordinates are outside the sheet matrix.
@test_throws ArgumentError fs(e + w_utm, n + h_utm)
fs(e + w_utm - 1, n + h_utm) == (2, 3)
@test_throws ArgumentError fs(e + w_utm, n + h_utm + 1)


##############################################################
# Test conversion from geographical coordinates to cell index
##############################################################


f1 = func_utm_to_cell_index(smb)
fi(e, n) = f1((e, n))
# The southwest corner, touched by four sheets, belongs to a sheet on the south side of the sheet matrix.
# Whose coordinates are outside the sheet matrix. We return the repeating index from the sheet, even
# if that one is outside the matrix
@test fi(e, n) == (1, 1)
# Our first sheet's top left
@test fi(e + ws_utm, n + hs_utm) == (1, 1)
# Our first sheets' bottom left
@test fi(e, n + 1)[1] == hs_utm / cell2utm
# Another utm coordinate with the same sheet index
@test fi(e, n + 2)[1] == hs_utm / cell2utm
# Our first sheets' bottom left
@test fi(e, n)[2] == 1
# Another utm coordinate with the same sheet index
@test fi(e + 1, n)[2] == 1
[fi(e, n + i)[1] for i in hs_utm:-1:1] == [1, 1, 2, 2, 3, 3, 4, 4]
[fi(e + i, n)[2] for i in 0:(ws_utm - 1)] == [1, 1, 2, 2, 3, 3]
