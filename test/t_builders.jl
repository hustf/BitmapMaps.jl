using Test
using BitmapMaps
using BitmapMaps: SheetMatrixBuilder, row_col_of_sheet

smb = SheetMatrixBuilder(3, 4, 1, 2, (43999, 6909048), 1, "nopath")
@test length(smb) == 6

(item, state) = iterate(smb) # 1, 2
(item, state) = iterate(smb, state) # 2, 3
(item, state) = iterate(smb, state) # 3, 4
(item, state) = iterate(smb, state) # 4, 5
(item, state) = iterate(smb, state) # 5, 6
(item, state) = iterate(smb, state) # 6, 7
next = iterate(smb, state)
@test isnothing(next)

# SheetMatrixBuilder corresponds to resource/matrix_sheet_pix_utm.svg

bm_pix_width, bm_pix_height, sheet_pix_width, sheet_pix_height, sw_corner, pix_to_utm_factor = (9, 8, 3, 4, (43999, 6909048), 2)

@test_throws ArgumentError SheetMatrixBuilder(bm_pix_width + 1, bm_pix_height, sheet_pix_width, sheet_pix_height, sw_corner, pix_to_utm_factor, "test")
@test_throws ArgumentError SheetMatrixBuilder(bm_pix_width, bm_pix_height + 1, sheet_pix_width, sheet_pix_height, sw_corner, pix_to_utm_factor, "test")
@test_throws ArgumentError SheetMatrixBuilder(bm_pix_width, bm_pix_height    , bm_pix_width * 2, sheet_pix_height, sw_corner, pix_to_utm_factor, "test")
@test_throws ArgumentError SheetMatrixBuilder(bm_pix_width, bm_pix_height, sheet_pix_width, bm_pix_height * 2, sw_corner, pix_to_utm_factor, "test")
smb = SheetMatrixBuilder(bm_pix_width, bm_pix_height, sheet_pix_width, sheet_pix_height, sw_corner, pix_to_utm_factor, "test")
@test length(smb) == 6
(item, state) = iterate(smb) # 1, 2
@test item.sheet_number == 1
pix_iter = item.pix_iter
@test length(pix_iter) == 12
@test pix_iter[1,1] == pix_iter[1]
@test pix_iter[2,1] == pix_iter[2]
@test pix_iter[3,1] == pix_iter[3]
f_I_to_utm = item.f_I_to_utm
@test_throws BoundsError f_I_to_utm(pix_iter[0, 1])
@test f_I_to_utm(pix_iter[1, 1]) .- sw_corner == (0, 8)
@test f_I_to_utm(pix_iter[2, 1]) .- sw_corner == (0, 6)
@test f_I_to_utm(pix_iter[3, 1]) .- sw_corner == (0, 4)
@test f_I_to_utm(pix_iter[4, 1]) .- sw_corner == (0, 2)
@test_throws BoundsError f_I_to_utm(pix_iter[5, 1])
@test_throws BoundsError f_I_to_utm(pix_iter[4, 0])
@test f_I_to_utm(pix_iter[4, 1]) .- sw_corner == (0, 2)
@test f_I_to_utm(pix_iter[4, 2]) .- sw_corner == (2, 2)
@test f_I_to_utm(pix_iter[4, 3]) .- sw_corner == (4, 2)
@test_throws BoundsError f_I_to_utm(pix_iter[4, 4])

@test state.sheet_number == 2
pix_iter = state.pix_iter
@test length(pix_iter) == 12
@test pix_iter[1,1] == pix_iter[1]
@test pix_iter[2,1] == pix_iter[2]
@test pix_iter[3,1] == pix_iter[3]
f_I_to_utm = state.f_I_to_utm
@test_throws BoundsError f_I_to_utm(pix_iter[0, 1])
@test f_I_to_utm(pix_iter[1, 1]) .- sw_corner == (0, 2 + 6 + 8)
@test f_I_to_utm(pix_iter[2, 1]) .- sw_corner == (0, 2 + 4 + 8)
@test f_I_to_utm(pix_iter[3, 1]) .- sw_corner == (0, 2 + 2 + 8)
@test f_I_to_utm(pix_iter[4, 1]) .- sw_corner == (0, 2 + 0 + 8)
@test_throws BoundsError f_I_to_utm(pix_iter[5, 1])
@test_throws BoundsError f_I_to_utm(pix_iter[4, 0])
@test f_I_to_utm(pix_iter[4, 1]) .- sw_corner == (0, 2  + 8)
@test f_I_to_utm(pix_iter[4, 2]) .- sw_corner == (2, 2 + 8)
@test f_I_to_utm(pix_iter[4, 3]) .- sw_corner == (4, 2 + 8)
@test_throws BoundsError f_I_to_utm(pix_iter[4, 4])


for (i, shp) in enumerate(smb)
     xy = shp.pixel_origin_ref_to_bitmapmap
     i == 1 && @test xy == (0, 4)
     i == 2 && @test xy == (0, 0)
     i == 3 && @test xy == (3, 4)
     i == 4 && @test xy == (3, 0)
     i == 5 && @test xy == (6, 4)
     i == 6 && @test xy == (6, 0)
     @test i == shp.sheet_number
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

for i in 1:6
    shp = smb[i]
    xy = shp.pixel_origin_ref_to_bitmapmap
    i == 1 && @test xy == (0, 4)
    i == 2 && @test xy == (0, 0)
    i == 3 && @test xy == (3, 4)
    i == 4 && @test xy == (3, 0)
    i == 5 && @test xy == (6, 4)
    i == 6 && @test xy == (6, 0)
    @test xy == BitmapMaps._SheetBuilder(smb, i).pixel_origin_ref_to_bitmapmap
end

for (i, shp) in enumerate(smb)
    shi1 = BitmapMaps._SheetBuilder(smb, i)
    @test shi1.f_I_to_utm(CartesianIndex(2, 2)) == shp.f_I_to_utm(CartesianIndex(2, 2))
end

@test_throws BoundsError smb[0,1]
@test_throws BoundsError smb[3,1]
@test_throws BoundsError smb[1,0]
@test_throws BoundsError smb[1,4]
@test smb[1,3] == smb[1, 3]
@test smb[1,3] == smb[5]
@test smb[2,3] == smb[6]
