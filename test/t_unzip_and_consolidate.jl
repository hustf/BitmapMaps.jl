using Test
using BitmapMaps
# SheetBuilder corresponds to resource/matrix_sheet_pix_utm.svg
bm_pix_width, bm_pix_height, sheet_pix_width, sheet_pix_height, sw_corner, pix_to_utm_factor = (9, 8, 3, 4, (43999, 6909048), 2)
smb = SheetMatrixBuilder(bm_pix_width, bm_pix_height, sheet_pix_width, sheet_pix_height, sw_corner, pix_to_utm_factor, "BitmapMaps/test")



for (p, zfi) in zip([smb[1,1], smb[2,2]], ["../resource/eksport_796345_20240420.zip", "../resource/eksport_796340_20240420.zip"])
    BitmapMaps.establish_folder(p)
    zipfi = joinpath(@__DIR__, zfi)
    dest = joinpath(BitmapMaps.full_folder_path(p), splitdir(zipfi)[2])
    if ! isfile(dest)
        cp(zipfi, dest)
    end
    # A downloaded and relevant zip file now exists as if downloaded by user.
    # Example export settings are found in /resource/hoydedata export settings.png
    # Other export options work, too.
    BitmapMaps.unzip_tif(p)
    if p == smb[2,2]
        @test isfile(joinpath(BitmapMaps.full_folder_path(p), "dom1", "data", "dom1-33-1-428-189-63.tif"))
    end
    if p == smb[1, 1]
        @test isfile(joinpath(BitmapMaps.full_folder_path(p), "dom1", "data", "dom1-33-1-428-189-53.tif"))
        @test isfile(joinpath(BitmapMaps.full_folder_path(p), "dom1", "data", "dom1-33-1-428-189-63.tif"))
    end
end




# Now consolidate (multiple) source files into CONSOLIDATED_FNAM in each sheet's folder.
for p in [smb[2,2], smb[1,1] ]    
    fnam_out = joinpath(BitmapMaps.full_folder_path(p), BitmapMaps.CONSOLIDATED_FNAM)
    if isfile(fnam_out)
        rm(fnam_out)
    end
    BitmapMaps.consolidate_elevation_data(p)
    @test isfile(fnam_out)
end

# Clean up
for p in [smb[2,2], smb[1,1] ]    
    fnam_out = joinpath(BitmapMaps.full_folder_path(p), BitmapMaps.CONSOLIDATED_FNAM)
    if isfile(fnam_out)
        rm(fnam_out)
    end
end


#=
# Do a bit of lower level testing....
import GeoArrays
using GeoArrays: crop, coords, bbox, Vertex, Center, indices
using GeoArrays: bbox!, GeoArray, AffineMap, SVector
import PNGFiles
using PNGFiles: Gray
#BitmapMaps.consolidate_elevation_data(smb)

p = smb[1,1]
fi = joinpath(BitmapMaps.full_folder_path(p), "dom1", "data", "dom1-33-1-428-189-63.tif")
r, c, min_x, min_y, max_x, max_y = BitmapMaps.parse_folder_name(BitmapMaps.full_folder_path(p))
cbox = (;min_x, min_y, max_x, max_y)
w = max_x - min_x
h = max_y - min_y

g = let
     A = zeros(Float32, w, h, 1)
     f = GeoArrays.AffineMap([1.0 0.0; 0.0 -1.0], 1.0 .* [min_x, max_y])
     GeoArray(A, f) 
end
# SW-NE Internal corner indices
sw_internal_indices = indices(g, SVector{2}(min_x, min_y + 1), Vertex()).I
ne_internal_indices = indices(g, SVector{2}(max_x - 1, max_y), Vertex()).I
@test sw_internal_indices == (1, 8)
@test ne_internal_indices == (6, 1)
min_i, min_j, max_i, max_j = BitmapMaps.indices_internal(g, min_x, min_y, max_x, max_y)


@test g.f(SVector(min_i - 1, min_j - 1)) .- (min_x, max_y) == [0, 0]
@test g.f(SVector(max_i    , min_j - 1)) .- (max_x, max_y) == [0, 0]
@test g.f(SVector(max_i    , min_j - 1)) .- (min_x, max_y) == [w, 0]
@test g.f(SVector(min_i - 1 + w, min_j - 1)) .- (min_x, max_y) == [w, 0]

@test integer_coords(g, min_i, min_j) == [min_x, max_y]
@test integer_coords(g, max_i, min_j) == [max_x - 1, max_y]
@test integer_coords(g, max_i, max_j) == [max_x - 1, min_y + 1]
@test integer_coords(g, min_i, max_j) == [min_x, min_y + 1]

bbox(g)






Int.(coords(g, CartesianIndex((min_i, min_j)), Vertex())) == (min_x, max_y)

coords(g, CartesianIndex((2, 2)), Vertex()) .- (min_x, min_y) .|> Int |> println
coords(g, CartesianIndex((1, 8)), Vertex()) .- (min_x, min_y) .|> Int |> println
coords(g, CartesianIndex((6, 8)), Vertex()) .- (min_x, min_y) .|> Int |> println











g.f(SVector(min_i, min_j))
coords(g,)





ii_min_x, ii_min_y = indices(g_dest, (cbox.min_x, cbox.min_y)).I
ii_max_x, ii_max_y = indices(g_dest, (cbox.max_x, cbox.max_y)).I
bbox!(g_dest, (;min_x = min_x -1, min_y = min_y - 1, max_x, max_y))
ii_min_x, ii_min_y = indices(g_dest, (cbox.min_x, cbox.min_y)).I
ii_max_x, ii_max_y = indices(g_dest, (cbox.max_x, cbox.max_y)).I
g_dest.f
g_dest.f.translation[2] - min_y
g_dest = let
    A = zeros(Float32, w, h, 1)
    f = GeoArrays.AffineMap([1.0 0.0; 0.0 -1.0], 1.0 .* [min_x , max_y ])
    GeoArray(A, f) 
end
@test length(collect(g_dest)) == w * h
@test cbox == bbox(g_dest)
@test first(eachindex(g_dest)) == CartesianIndex(1, 1, 1)
@test last(eachindex(g_dest)) == CartesianIndex(w, h, 1)

@test indices(g_dest, (min_x, max_y), Vertex()).I == (1, 1)
indices(g_dest, (min_x, max_y - 1), Vertex()).I
indices(g_dest, SVector{2}(min_x, min_y + 1), Vertex()).I





coords(g_dest, CartesianIndex((2, 2)), Vertex()) .- (min_x, min_y) .|> Int |> println
coords(g_dest, CartesianIndex((1, 8)), Vertex()) .- (min_x, min_y) .|> Int |> println
coords(g_dest, CartesianIndex((6, 8)), Vertex()) .- (min_x, min_y) .|> Int |> println


ii_min_x, ii_min_y = indices(g_dest, (cbox.min_x, cbox.min_y), Vertex()).I
ii_max_x, ii_max_y = indices(g_dest, (cbox.max_x, cbox.max_y), Vertex()).I

for ii in eachindex(g_dest)
    @show ii.I
end





ga = GeoArrays.read(fi)
#
indices(ga, southwest_corner(p), Vertex())
indices(ga, (min_x, max_y), Vertex())  
indices(ga, (min_x + 1, min_y), Vertex())  
coords(ga, (1, 553), Vertex()) == [min_x, min_y]
@test ga[indices(ga, (min_x + 1, max_y), Vertex())] > 0.0
# It unfortunately seems this file haven't got data for min_x. 
fi = joinpath(BitmapMaps.full_folder_path(p), "dom1", "data", "dom1-33-1-428-189-53.tif")
ga = GeoArrays.read(fi)
indices(ga, southwest_corner(p), Vertex())
indices(ga, (min_x, max_y), Vertex())  
indices(ga, (min_x + 1, min_y), Vertex())  
coords(ga, (800, 553), Vertex()) == [min_x, min_y]
# This file contains min_x
@test ga[indices(ga, (min_x, max_y), Vertex())] > 0.0

coords(ga, (1, 1), Vertex())






indices(ga, (min_x + 2, min_y), Vertex())  
indices(ga, (min_x + 3, min_y), Vertex())  
indices(ga, (min_x + 4, min_y), Vertex())  
indices(ga, (min_x + 5, min_y), Vertex())  
indices(ga, (min_x + w, min_y), Vertex())  
indices(ga, (max_x, max_y), Vertex())
indices(ga, (min_x + w + 1, min_y), Vertex())  
@test 
indices(ga, southwest_corner(p), Center())
indices(ga, (min_x, max_y), Center())  
indices(ga, (min_x + 1, max_y), Center())  
indices(ga, (min_x + 2, max_y), Center())  
indices(ga, (min_x + 3, max_y), Center())  
indices(ga, (min_x + 4, max_y), Center())  
indices(ga, (min_x + 5, max_y), Center())  
indices(ga, (min_x + w, max_y), Center())  
indices(ga, (min_x + w + 1, max_y), Center())  
indices(ga, (max_x, max_y), Center())  

indices(ga, southwest_corner(p))
indices(ga, (min_x, max_y))  
indices(ga, (min_x + 1, max_y))  
indices(ga, (min_x + 2, max_y))  
indices(ga, (min_x + 3, max_y))  
indices(ga, (min_x + 4, max_y))  
indices(ga, (min_x + 5, max_y))  
indices(ga, (min_x + w, max_y))  
indices(ga, (min_x + w + 1, max_y), Vertex())  
indices(ga, (max_x, max_y))  

# Conclusion: Always use Vertex(). Don't use crop at all.


indices(ga, northeast_external_corner(p), Vertex())
indices(ga, (max_x, min_y), Vertex())
ga[1, 553]

@test BitmapMaps.boxwidth(g_source) == w
# The crop function returned a taller box than wanted. 
BitmapMaps.boxheight(g_source) == h + 1
# ...but which row should we ignore?
display(map(pix -> Gray(max(0.0, pix) / 1500), g_source.A[:,:, 1]))
display(map(pix -> Gray(max(0.0, pix) / 1500), ga.A[:,:, 1]))
g_source.A[0, 0]

Tuple(indices(g_source, southwest_corner(p), Vertex())) 
Tuple(indices(g_source, southwest_corner(p), Center())) 
Tuple(indices(ga, southwest_corner(p), Vertex())) 
Tuple(indices(ga, southwest_corner(p), Center())) 

Tuple(indices(g_source, northeast_internal_corner(p), Vertex())) 
Tuple(indices(g_source, northeast_external_corner(p), Vertex())) 
Tuple(indices(ga, northeast_external_corner(p), Vertex())) 



coords(g_source)[1]



Tuple(indices(g_source, nw_source, Vertex())) # TODO add Vertex
=#