using Test
using BitmapMaps
using BitmapMaps: BmPartition

bmp = BmPartition(3, 4, 1, 2, (43999, 6909048), 1)
@test length(bmp) == 6

(item, state) = iterate(bmp) # 1, 2
(item, state) = iterate(bmp, state) # 2, 3
(item, state) = iterate(bmp, state) # 3, 4
(item, state) = iterate(bmp, state) # 4, 5
(item, state) = iterate(bmp, state) # 5, 6
(item, state) = iterate(bmp, state) # 6, 7
next = iterate(bmp, state)
@test isnothing(next)

# Corresponding to resource/map_sheet_utm_pix.svg

bm_pixel_width, bm_pixel_height, sheet_width, sheet_height, bm_southwest_corner, pixel_distance = (9, 8, 3, 4, (43999, 6909048), 2)
bmp = BmPartition(bm_pixel_width, bm_pixel_height, sheet_width, sheet_height, bm_southwest_corner, pixel_distance)
@test length(bmp) == 6
(item, state) = iterate(bmp) # 1, 2
@test item[2] == 1
pix_iter = item[1].pix_iter
@test length(pix_iter) == 12
@test pix_iter[1,1] == pix_iter[1]
@test pix_iter[2,1] == pix_iter[2]
@test pix_iter[3,1] == pix_iter[3]
f_I_to_utm = item[1].f_I_to_utm
@test_throws BoundsError f_I_to_utm(pix_iter[0, 1])
@test f_I_to_utm(pix_iter[1, 1]) .- bm_southwest_corner == (0, 6)
@test f_I_to_utm(pix_iter[2, 1]) .- bm_southwest_corner == (0, 4)
@test f_I_to_utm(pix_iter[3, 1]) .- bm_southwest_corner == (0, 2)
@test f_I_to_utm(pix_iter[4, 1]) .- bm_southwest_corner == (0, 0)
@test_throws BoundsError f_I_to_utm(pix_iter[5, 1])
@test_throws BoundsError f_I_to_utm(pix_iter[4, 0])
@test f_I_to_utm(pix_iter[4, 1]) .- bm_southwest_corner == (0, 0)
@test f_I_to_utm(pix_iter[4, 2]) .- bm_southwest_corner == (2, 0)
@test f_I_to_utm(pix_iter[4, 3]) .- bm_southwest_corner == (4, 0)
@test_throws BoundsError f_I_to_utm(pix_iter[4, 4])

@test state[2] == 2
pix_iter = state[1].pix_iter
@test length(pix_iter) == 12
@test pix_iter[1,1] == pix_iter[1]
@test pix_iter[2,1] == pix_iter[2]
@test pix_iter[3,1] == pix_iter[3]
f_I_to_utm = state[1].f_I_to_utm
@test_throws BoundsError f_I_to_utm(pix_iter[0, 1])
@test f_I_to_utm(pix_iter[1, 1]) .- bm_southwest_corner == (0, 6 + 8)
@test f_I_to_utm(pix_iter[2, 1]) .- bm_southwest_corner == (0, 4 + 8)
@test f_I_to_utm(pix_iter[3, 1]) .- bm_southwest_corner == (0, 2 + 8)
@test f_I_to_utm(pix_iter[4, 1]) .- bm_southwest_corner == (0, 0 + 8)
@test_throws BoundsError f_I_to_utm(pix_iter[5, 1])
@test_throws BoundsError f_I_to_utm(pix_iter[4, 0])
@test f_I_to_utm(pix_iter[4, 1]) .- bm_southwest_corner == (0, 0  + 8)
@test f_I_to_utm(pix_iter[4, 2]) .- bm_southwest_corner == (2, 0 + 8)
@test f_I_to_utm(pix_iter[4, 3]) .- bm_southwest_corner == (4, 0 + 8)
@test_throws BoundsError f_I_to_utm(pix_iter[4, 4])


for (i, shp) in enumerate(bmp)
     display(shp[1])
     xy = shp[1].pixel_origin_ref_to_bitmapmap
     i == 1 && @test xy == (0, 4)
     i == 2 && @test xy == (0, 0)
     i == 3 && @test xy == (3, 4)
     i == 4 && @test xy == (3, 0)
     i == 5 && @test xy == (6, 4)
     i == 6 && @test xy == (6, 0)
     @test i == shp[2]
end

for (i, _) in enumerate(bmp)
    current_row, current_col = row_col_of_sheet(bmp, i)
    i == 1 && @test (current_row, current_col) == (1, 1)
    i == 2 && @test (current_row, current_col) == (2, 1)
    i == 3 && @test (current_row, current_col) == (1, 2)
    i == 4 && @test (current_row, current_col) == (2, 2)
    i == 5 && @test (current_row, current_col) == (1, 3)
    i == 6 && @test (current_row, current_col) == (2, 3)
end

for i in 1:6
    shp = bmp[i]
    xy = shp[1].pixel_origin_ref_to_bitmapmap
    i == 1 && @test xy == (0, 4)
    i == 2 && @test xy == (0, 0)
    i == 3 && @test xy == (3, 4)
    i == 4 && @test xy == (3, 0)
    i == 5 && @test xy == (6, 4)
    i == 6 && @test xy == (6, 0)
    @test xy == SheetPartition(bmp, i).pixel_origin_ref_to_bitmapmap
end

for (i, shp) in enumerate(bmp)
    shi1 = SheetPartition(bmp, i)
    @show i
    display(shp[1])
    display(shi1)
    @test shi1.f_I_to_utm(CartesianIndex(2, 2)) == shp[1].f_I_to_utm(CartesianIndex(2, 2))
end



SheetPartition(bmp, 1)
SheetPartition(bmp, 2)
SheetPartition(bmp, 3)


for (i, shp) in enumerate(bmp)
    display(shp[1])
end