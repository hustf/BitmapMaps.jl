using Test
using BitmapMaps
using BitmapMaps: BmPartition

bmp = BmPartition(3, 4, 1, 2, (43999, 6909048), 1, "nopath")
@test length(bmp) == 6

(item, state) = iterate(bmp) # 1, 2
(item, state) = iterate(bmp, state) # 2, 3
(item, state) = iterate(bmp, state) # 3, 4
(item, state) = iterate(bmp, state) # 4, 5
(item, state) = iterate(bmp, state) # 5, 6
(item, state) = iterate(bmp, state) # 6, 7
next = iterate(bmp, state)
@test isnothing(next)

# BM paritioning corresponds to resource/map_sheet_utm_pix.svg

pix_width, pix_height, sheet_pix_width, sheet_pix_height, sw_corner, pix_to_utm_factor = (9, 8, 3, 4, (43999, 6909048), 2)
bmp = BmPartition(pix_width, pix_height, sheet_pix_width, sheet_pix_height, sw_corner, pix_to_utm_factor, "test")
@test length(bmp) == 6
(item, state) = iterate(bmp) # 1, 2
@test item.sheet_number == 1
pix_iter = item.pix_iter
@test length(pix_iter) == 12
@test pix_iter[1,1] == pix_iter[1]
@test pix_iter[2,1] == pix_iter[2]
@test pix_iter[3,1] == pix_iter[3]
f_I_to_utm = item.f_I_to_utm
@test_throws BoundsError f_I_to_utm(pix_iter[0, 1])
@test f_I_to_utm(pix_iter[1, 1]) .- sw_corner == (0, 6)
@test f_I_to_utm(pix_iter[2, 1]) .- sw_corner == (0, 4)
@test f_I_to_utm(pix_iter[3, 1]) .- sw_corner == (0, 2)
@test f_I_to_utm(pix_iter[4, 1]) .- sw_corner == (0, 0)
@test_throws BoundsError f_I_to_utm(pix_iter[5, 1])
@test_throws BoundsError f_I_to_utm(pix_iter[4, 0])
@test f_I_to_utm(pix_iter[4, 1]) .- sw_corner == (0, 0)
@test f_I_to_utm(pix_iter[4, 2]) .- sw_corner == (2, 0)
@test f_I_to_utm(pix_iter[4, 3]) .- sw_corner == (4, 0)
@test_throws BoundsError f_I_to_utm(pix_iter[4, 4])

@test state.sheet_number == 2
pix_iter = state.pix_iter
@test length(pix_iter) == 12
@test pix_iter[1,1] == pix_iter[1]
@test pix_iter[2,1] == pix_iter[2]
@test pix_iter[3,1] == pix_iter[3]
f_I_to_utm = state.f_I_to_utm
@test_throws BoundsError f_I_to_utm(pix_iter[0, 1])
@test f_I_to_utm(pix_iter[1, 1]) .- sw_corner == (0, 6 + 8)
@test f_I_to_utm(pix_iter[2, 1]) .- sw_corner == (0, 4 + 8)
@test f_I_to_utm(pix_iter[3, 1]) .- sw_corner == (0, 2 + 8)
@test f_I_to_utm(pix_iter[4, 1]) .- sw_corner == (0, 0 + 8)
@test_throws BoundsError f_I_to_utm(pix_iter[5, 1])
@test_throws BoundsError f_I_to_utm(pix_iter[4, 0])
@test f_I_to_utm(pix_iter[4, 1]) .- sw_corner == (0, 0  + 8)
@test f_I_to_utm(pix_iter[4, 2]) .- sw_corner == (2, 0 + 8)
@test f_I_to_utm(pix_iter[4, 3]) .- sw_corner == (4, 0 + 8)
@test_throws BoundsError f_I_to_utm(pix_iter[4, 4])


for (i, shp) in enumerate(bmp)
     display(shp)
     xy = shp.pixel_origin_ref_to_bitmapmap
     i == 1 && @test xy == (0, 4)
     i == 2 && @test xy == (0, 0)
     i == 3 && @test xy == (3, 4)
     i == 4 && @test xy == (3, 0)
     i == 5 && @test xy == (6, 4)
     i == 6 && @test xy == (6, 0)
     @test i == shp.sheet_number
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
    xy = shp.pixel_origin_ref_to_bitmapmap
    i == 1 && @test xy == (0, 4)
    i == 2 && @test xy == (0, 0)
    i == 3 && @test xy == (3, 4)
    i == 4 && @test xy == (3, 0)
    i == 5 && @test xy == (6, 4)
    i == 6 && @test xy == (6, 0)
    @test xy == BitmapMaps._SheetPartition(bmp, i).pixel_origin_ref_to_bitmapmap
end

for (i, shp) in enumerate(bmp)
    shi1 = BitmapMaps._SheetPartition(bmp, i)
    @show i
    display(shp)
    display(shi1)
    @test shi1.f_I_to_utm(CartesianIndex(2, 2)) == shp.f_I_to_utm(CartesianIndex(2, 2))
end

@test_throws BoundsError bmp[0,1]
@test_throws BoundsError bmp[3,1]
@test_throws BoundsError bmp[1,0]
@test_throws BoundsError bmp[1,4]
@test bmp[1,3] == bmp[1, 3]
@test bmp[1,3] == bmp[5]
@test bmp[2,3] == bmp[6]

@test southwest_corner(bmp) == southwest_corner(bmp[1, 1]) 
@test northeast_external_corner(bmp) == northeast_external_corner(bmp[2, 3])
@test northeast_internal_corner(bmp) == northeast_internal_corner(bmp[2, 3])

@test geo_centre(bmp) isa Tuple{Float64, Float64}
@test geo_centre(bmp[1,1]) isa Tuple{Float64, Float64}
@test geo_grid_centre_single(bmp) isa Tuple{Int64, Int64}
@test geo_grid_centre_single(bmp[1,1]) isa Tuple{Int64, Int64}


@test all(geo_centre(bmp)                  .> southwest_corner(bmp))
@test all(geo_centre(bmp[1,1])             .> southwest_corner(bmp[1,1]))
@test all(geo_grid_centre_single(bmp)      .> southwest_corner(bmp))
@test all(geo_grid_centre_single(bmp[1,1]) .> southwest_corner(bmp[1,1]))

@test all(geo_centre(bmp)                  .< northeast_internal_corner(bmp))
@test all(geo_centre(bmp[1,1])             .< northeast_internal_corner(bmp[1,1]))
@test all(geo_grid_centre_single(bmp)      .< northeast_internal_corner(bmp))
@test all(geo_grid_centre_single(bmp[1,1]) .< northeast_internal_corner(bmp[1,1]))


# Check that a sheet's folder name matches its northeast_external_corner
@test tryparse.(Int, split(bmp[2,2].pthsh, ' ')[end - 1: end]) == collect(northeast_external_corner(bmp[2, 2]))


@test geo_area(bmp) == bmp.pix_width * bmp. pix_height * pix_to_utm_factor^2
@test geo_area(bmp[1]) == (bmp.pix_width * bmp. pix_height * pix_to_utm_factor^2) / (bmp.nrows * bmp.ncols)

@test bounding_box_external_string(bmp) == "(43999 6909048)-(44017 6909064)"
@test bounding_box_external_string(bmp[1,1]) == "(43999 6909048)-(44005 6909056)"
@test replace(BitmapMaps.bounding_box_polygon_string(bmp), "\t" => "") == """
    POLYGON ((43999 6909048, 44005 6909048, 44005 6909056, 43999 6909056, 43999 6909048),
    (43999 6909056, 44005 6909056, 44005 6909064, 43999 6909064, 43999 6909056),
    (44005 6909048, 44011 6909048, 44011 6909056, 44005 6909056, 44005 6909048),
    (44005 6909056, 44011 6909056, 44011 6909064, 44005 6909064, 44005 6909056),
    (44011 6909048, 44017 6909048, 44017 6909056, 44011 6909056, 44011 6909048),
    (44011 6909056, 44017 6909056, 44017 6909064, 44011 6909064, 44011 6909056))"""