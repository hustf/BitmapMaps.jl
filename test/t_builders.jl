using Test
using BitmapMaps
using BitmapMaps: SheetMatrixBuilder, row_col_of_sheet


smb = SheetMatrixBuilder((43999, 6909048), (2, 3), 1, 191, 275, 11811, "nopath")
@test length(smb) == 6

(item, state) = iterate(smb) # 1, 2
(item, state) = iterate(smb, state) # 2, 3
(item, state) = iterate(smb, state) # 3, 4
(item, state) = iterate(smb, state) # 4, 5
(item, state) = iterate(smb, state) # 5, 6
(item, state) = iterate(smb, state) # 6, 7
next = iterate(smb, state)
@test isnothing(next)


# The alternative constructor
# Corresponds to resource/matrix_sheet_cell_utm.svg
sw_corner, nrc, cell2utm = (44000, 6909047), (2,3), 2
sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, pth = 3, 4, 1000, "nopath"

@test_throws MethodError SheetMatrixBuilder(sw_corner .* 1.0, nrc, cell2utm, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, pth)
@test_throws MethodError SheetMatrixBuilder(sw_corner, nrc .* 1.0, cell2utm, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, pth)
@test_throws "cell_to_utm_factor" SheetMatrixBuilder(sw_corner, nrc, cell2utm - 2, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, pth)
@test_throws "sheet_width_mm" SheetMatrixBuilder(sw_corner, nrc, cell2utm, -sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, pth)
@test_throws "sheet_height_mm" SheetMatrixBuilder(sw_corner, nrc, cell2utm, sheet_width_mm, -sheet_height_mm, density_pt_m⁻¹, pth)
@test_throws "density_pt_m⁻¹" SheetMatrixBuilder(sw_corner, nrc, cell2utm, sheet_width_mm, sheet_height_mm, -density_pt_m⁻¹, pth)
@test_throws "pth" SheetMatrixBuilder(sw_corner, nrc, cell2utm, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, "")
smb = SheetMatrixBuilder(sw_corner, nrc, cell2utm, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, pth)

@test length(smb) == 6
(item, state) = iterate(smb) # 1, 2
@test item.sheet_number == 1
cell_iter = item.cell_iter
@test length(cell_iter) == 12
@test cell_iter[1,1] == cell_iter[1]
@test cell_iter[2,1] == cell_iter[2]
@test cell_iter[3,1] == cell_iter[3]
f_I_to_utm = item.f_I_to_utm
@test_throws BoundsError f_I_to_utm(cell_iter[0, 1])
@test f_I_to_utm(cell_iter[1, 1]) .- sw_corner == (0, 8)
@test f_I_to_utm(cell_iter[2, 1]) .- sw_corner == (0, 6)
@test f_I_to_utm(cell_iter[3, 1]) .- sw_corner == (0, 4)
@test f_I_to_utm(cell_iter[4, 1]) .- sw_corner == (0, 2)
@test_throws BoundsError f_I_to_utm(cell_iter[5, 1])
@test_throws BoundsError f_I_to_utm(cell_iter[4, 0])
@test f_I_to_utm(cell_iter[4, 1]) .- sw_corner == (0, 2)
@test f_I_to_utm(cell_iter[4, 2]) .- sw_corner == (2, 2)
@test f_I_to_utm(cell_iter[4, 3]) .- sw_corner == (4, 2)
@test_throws BoundsError f_I_to_utm(cell_iter[4, 4])

@test state.sheet_number == 2
cell_iter = state.cell_iter
@test length(cell_iter) == 12
@test cell_iter[1,1] == cell_iter[1]
@test cell_iter[2,1] == cell_iter[2]
@test cell_iter[3,1] == cell_iter[3]
f_I_to_utm = state.f_I_to_utm
@test_throws BoundsError f_I_to_utm(cell_iter[0, 1])
@test f_I_to_utm(cell_iter[1, 1]) .- sw_corner == (0, 2 + 6 + 8)
@test f_I_to_utm(cell_iter[2, 1]) .- sw_corner == (0, 2 + 4 + 8)
@test f_I_to_utm(cell_iter[3, 1]) .- sw_corner == (0, 2 + 2 + 8)
@test f_I_to_utm(cell_iter[4, 1]) .- sw_corner == (0, 2 + 0 + 8)
@test_throws BoundsError f_I_to_utm(cell_iter[5, 1])
@test_throws BoundsError f_I_to_utm(cell_iter[4, 0])
@test f_I_to_utm(cell_iter[4, 1]) .- sw_corner == (0, 2  + 8)
@test f_I_to_utm(cell_iter[4, 2]) .- sw_corner == (2, 2 + 8)
@test f_I_to_utm(cell_iter[4, 3]) .- sw_corner == (4, 2 + 8)
@test_throws BoundsError f_I_to_utm(cell_iter[4, 4])


for (i, sb) in enumerate(smb)
     xy = sb.pixel_origin_ref_to_bitmapmap
     i == 1 && @test xy == (0, 4)
     i == 2 && @test xy == (0, 0)
     i == 3 && @test xy == (3, 4)
     i == 4 && @test xy == (3, 0)
     i == 5 && @test xy == (6, 4)
     i == 6 && @test xy == (6, 0)
     @test i == sb.sheet_number
end

for (i, _) in enumerate(smb)
    current_row, current_col = row_col_of_sheet(smb, i)
    i == 1 && @test (current_row, current_col) == (1, 1)
    i == 2 && @test (current_row, current_col) == (2, 1)
    i == 3 && @test (current_row, current_col) == (1, 2)
    i == 4 && @test (current_row, current_col) == (2, 2)
    i == 5 && @test (current_row, current_col) == (1, 3)
    i == 6 && @test (current_row, current_col) == (2, 3)
end

let
    for i in 1:6
        sb = smb[i]
        xy = sb.pixel_origin_ref_to_bitmapmap
        i == 1 && @test xy == (0, 4)
        i == 2 && @test xy == (0, 0)
        i == 3 && @test xy == (3, 4)
        i == 4 && @test xy == (3, 0)
        i == 5 && @test xy == (6, 4)
        i == 6 && @test xy == (6, 0)
        @test xy == BitmapMaps._SheetBuilder(smb, i).pixel_origin_ref_to_bitmapmap
    end
end
for (i, sb) in enumerate(smb)
    shi1 = BitmapMaps._SheetBuilder(smb, i)
    @test shi1.f_I_to_utm(CartesianIndex(2, 2)) == sb.f_I_to_utm(CartesianIndex(2, 2))
end

@test_throws BoundsError smb[0,1]
@test_throws BoundsError smb[3,1]
@test_throws BoundsError smb[1,0]
@test_throws BoundsError smb[1,4]
@test smb[1,3] == smb[1, 3]
@test smb[1,3] == smb[5]
@test smb[2,3] == smb[6]

@test repr(smb[2,3]) == "SheetBuilder((6, 0), (1:4, 1:3), (44012, 6909063)@(1, 1), 6, 1000, \"nopath\\\\2-3__44012-6909055__44018-6909063\")\n"
@test repr(smb) == "SheetMatrixBuilder((44000, 6909047), CartesianIndices((1:2, 1:3)), 2, 3, 4, 1000, \"nopath\")"

