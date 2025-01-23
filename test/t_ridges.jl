# This creates some ridges and dieders which would cross all four boundaries
# Visually check that we don't draw lines parallell to the border 
# on the bord

using Test
using BitmapMaps
using BitmapMaps: smooth_laplacian, _corners_and_dieders!, mark_at!, line!, indices_of_border
if ! @isdefined mountain
    include("t_contour_func.jl")
end
side = 1000
si = CartesianIndices((1:1:side, 1:1:side))
z = Matrix{Gray{Float32}}(undef, side, side)
fill!(z, -10f0)
for y = 1:100:500
    for l in 1:1
        line!(z, CartesianIndex(1, 3y + l), CartesianIndex(y + l, 1),  √2 )
        line!(z, CartesianIndex(1000, 1000 - y - l), CartesianIndex(1000 - 2y - l, 1000),  √2 )
    end
end
z .+= 20.0f0
z +=  test_terrain(;P = Matrix{Gray{Float32}}(undef, side, side))
display_if_vscode(z)
# Prepare buffers
result = zeros(RGB{N0f8}, side, side)
fill!(result, RGB{N0f8}(0.5, 0.4, 0.1))
bbuf = Array{Gray{Bool}}(undef, side, side)
fill!(bbuf, false)
masked_laplacian = smooth_laplacian(z, si) #.* (1 .- bumpy_patch(z, si))
extrema(masked_laplacian)
display_if_vscode(masked_laplacian)
colo_corner = RGBA{N0f8}(0.118, 0.102, 0.141, 0.85)
colo_dieder = RGBA{N0f8}(0.0, 0.64, 1.0, 0.7)  
criterion_functions = [<(-0.02f0), >(0.05f0)]
thicknesses = [3, 5]
_corners_and_dieders!(result, bbuf, masked_laplacian, criterion_functions, [colo_corner, colo_dieder], thicknesses)
display_if_vscode(result)
