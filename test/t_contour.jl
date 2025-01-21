# This optimises parameter and algorithms for identifying forest areas.
#
# This is best done in a slightly large area, which we won't 
# include in the repo resources. Instead, we rely on local files,
# which can be downloaded without creating a user account.
using Test
using BitmapMaps
using BitmapMaps: fir_lp_coefficients, centered
using BitmapMaps: imfilter, mapwindow, mapwindow!, define_builder
using BitmapMaps: _elev_contours, __elev_contours, Gray, load, felzenszwalb
using BitmapMaps: contour_lines_overlay, CONTOUR_FNAM, CONSOLIDATED_FNAM, TOPORELIEF_FNAM, COMPOSITE_FNAM
using BitmapMaps: RGB, RGBA, labels_map, N0f8, segment_mean, segment_pixel_count, dilate!
using BitmapMaps: CompositeDestinationOver, cartesian_index_string, imfilter, bumpy_patch
import BitmapMaps: bumpy_patch, bumpiness, blue, red, green, smooth_surface_fir, roughness_of_surface_fir
using BitmapMaps: cell_to_utm_factor
import ImageSegmentation
using ImageSegmentation: unseeded_region_growing, region_splitting, fast_scanning, prune_segments
using Statistics: quantile
using Random
#=
include("t_contour_func.jl")


smb = define_builder(;pth = "BitmapMaps/1 1 18145 6894028 71625 6971028", 
    cell_to_utm_factor = 1, 
    southwest_corner = (35683, 6919318),
    density_pt_m⁻¹ = 11811,
    nrc = (8, 10)
    )
for i in 1:6
    @test ispath(full_folder_path(smb[i])) # If not:  run_bitmapmap_pipeline(smb)
end









#
# A test image
#
 
z =  test_terrain()

# topo relief, for overview
topo = RGB{N0f8}.(BitmapMaps.__topo_relief(Float32.(z), CartesianIndices(size(z)), 1))
display_if_vscode(topo[900:1050, 355:1065])

# Smoothing, to be used in foresty areas:
Δz_ideal = .- forest(z) .- noise(z)
Δz_smooth = z .- smooth_surface_fir(z; w = 49)
-Δz_smooth
-Δz_smooth[900:1050, 355:1065]
-Δz_ideal[900:1050, 355:1065]
#
# Optimize  smooth_surface_fir parameter w
#
f_evaluate_fir = w -> let 
    Δz_smooth = z .- smooth_surface_fir(z; w)
    Float64(sum(abs.(Δz_smooth .- Δz_ideal)))
end
find_minimimzing_param(f_evaluate_fir, 27:2:199)
# Check minimum
f_evaluate_fir(101) - f_evaluate_fir(99)
f_evaluate_fir(103) - f_evaluate_fir(101)






# Optimize parameters for roughness_of_surface_fir
#

for window in 41:2:85
    z_rough = roughness_of_surface_fir(z; w = window)
    @show window
    display(z_rough[900:1050, 355:1065])
    if @isdefined UnicodePlots
        p = lineplot(Float32.(Δz_ideal[600:1050, 900]); width = 200, height = 30, title = string(window), ylim = (-8, 3))
        p = lineplot!(p, Float32.(z_rough[600:1050, 900]))
        println(p)
    end
    sleep(0.5)
end


for window in 49:2:49
    z_bump = bumpiness(z; w = window)
    display(z_bump[900:1050, 355:1065])
    sleep(1)
    z_bump1 = bumpy_patch(z; w = 49, amp_cut⁻ = 1.0, amp_cut⁺ =  4.6477094f0, z_max = 680.0f0)
    display(z_bump1[800:1050, 355:1065])
    @show window
    if @isdefined UnicodePlots
        p = lineplot(Float32.(Δz_ideal[600:1050, 800]); width = 200, height = 30, title = string(window), ylim = (-8, 3))
        p = lineplot!(p, Float32.(z_bump[600:1050, 800]))
        p = lineplot!(p, Float32.(z_bump1[600:1050, 800]))
        println(p)
    end
    sleep(0.5)
end



vthick = [1, 3, 5]
vdist = [20, 100, 1000] # Contour spacings.
minlen = 6 
sb = smb[3, 3]

g = readclose(joinpath(full_folder_path(sb), CONSOLIDATED_FNAM))
ny, nx = size(sb.cell_iter)
source_indices = (1:cell_to_utm_factor(sb):(ny  * cell_to_utm_factor(sb)), 1:cell_to_utm_factor(sb):(nx * cell_to_utm_factor(sb)))
si = CartesianIndices(source_indices)
sicr = CartesianIndices((4500:1:(4500  * cell_to_utm_factor(sb)), 2000:1:(nx * cell_to_utm_factor(sb))))
bump = bumpiness(transpose(g.A[:, :, 1])[sicr])
inds = CartesianIndices((1750:2050, 350:1100))
bump[inds]
map(x -> Gray{N0f8}(x > 0.6), smooth_surface_fir(bump; nyquist_denom = 4, w = 49))[inds]
map(x -> Gray{N0f8}(x > 0.5), smooth_surface_fir(bump; nyquist_denom = 4, w = 49))[inds]
map(x -> Gray{N0f8}(x > 0.4), smooth_surface_fir(bump; nyquist_denom = 4, w = 49))[inds]
map(x -> Gray{N0f8}(x > 0.3), smooth_surface_fir(bump; nyquist_denom = 4, w = 49))[inds]
map(x -> Gray{N0f8}(x > 0.2), smooth_surface_fir(bump; nyquist_denom = 4, w = 49))[inds]
map(x -> Gray{N0f8}(x > 0.1), smooth_surface_fir(bump; nyquist_denom = 4, w = 49))[inds]
map(x -> Gray{N0f8}(x > 0.1), smooth_surface_fir(bump; nyquist_denom = 4, w = 49))[inds]

d = map(x -> Gray{N0f8}(x > 0.6), smooth_surface_fir(bump; nyquist_denom = 4, w = 49))
d[inds]

cell_count_min = 9062
for k in 0.2:0.1:1.0
    gg = fast_scanning(d, k)
    f_to_be_removed = i -> (segment_pixel_count(gg, i) < cell_count_min)
    f_diff(rem_label, neigh_label) = segment_mean(gg, rem_label) - segment_mean(gg, neigh_label)
    h = prune_segments(gg, f_to_be_removed, f_diff)
    img = map(i -> get_random_color(i), labels_map(h))
    open_as_temp_in_gimp(img)
    @show k length(gg.segment_labels) length(h.segment_labels)
end




begin
    nyquist_denom = 4
    inds = CartesianIndices((1750:1750, 350:1100))
    p = lineplot(vec((Float32.(bump[inds]))); width = 200, height = 25)
    d = smooth_surface_fir(bump; nyquist_denom, w = 15)[inds]
    @show count(x-> x> 0, vec(Float32.(d)))
    p = lineplot!(p, vec(Float32.(d)))
    #p = lineplot(Float32.(d[600:1050, 1500]))
end

mask = bumpy_patch(transpose(g.A[:, :, 1])[sicr])
open_as_temp_in_gimp(mask)
mask
conto = _elev_contours(full_folder_path(sb), sb.cell_iter, cell_to_utm_factor(sb), minlen, vthick, vdist)
open_as_temp_in_gimp(conto)

















# Pre-process as in _elev_contours
bh = bumpy_patch(z)
@show extrema(bh)
display_if_vscode(bh[800:1050, 355:1065])
z_smooth = smooth_surface_fir(z)
display_if_vscode(z_smooth[800:1050, 355:1065])

img = RGB{Float32}.(bh, # red
        z,           # green
        Float32.(z_smooth)) # blue
# Test contours
conto = __elev_contours(img, minlen, vthick, vdist)


fnam = "tempo.png"
save_png_with_phys(fnam, conto)
gimp_path = "C:\\Program Files\\GIMP 2\\bin\\gimp-2.10.exe"
@async run(`$gimp_path "$fnam"`)


fnam = "temp.png"
save_png_with_phys(fnam, topo)
gimp_path = "C:\\Program Files\\GIMP 2\\bin\\gimp-2.10.exe"
@async run(`$gimp_path "$fnam"`)

sb = smb[2,2]
img = _elev_contours(full_folder_path(sb), sb.cell_iter, cell_to_utm_factor(sb), minlen, vthick, vdist)














# DELETE

# The low-pass filter deforms the large-scale terrain a lot. We must do better.
topo = RGB{N0f8}.(BitmapMaps.__topo_relief(Float32.(z), CartesianIndices(size(z)), 1))
display_if_vscode(topo[600:1050, 355:765])



for window in 41:2:85
    z_rough = roughness_of_surface_fir(z; w = window)
    @show window
    display(z_rough[900:1050, 355:1065])
    p = lineplot(Float32.(Δz_ideal[600:1050, 900]); width = 200, height = 30, title = string(window), ylim = (-8, 3))
    p = lineplot!(p, Float32.(z_rough[600:1050, 900]))
    println(p)
    sleep(0.5)
 end


for window in 3:2:49
   Δz_smooth = z .- smooth_surface_fir(z; w = window)
   @show window
   display_if_vscode(-Δz_smooth[900:1050, 355:1065])
   sleep(0.5)
end


    @show window extrema(Float32.(Δz))
    p = lineplot(Float32.(z[350, 800:1200]); width = 200, height = 50)
    p = lineplot!(p, Float32.(z_smooth[350, 800:1200]))
    println(p)
    sleep(1)


    devi = Δz .- Δz_ideal
    @show window extrema(Float32.(devi))
    -devi[300:400, 55:1850]
    -devi[300:400, 600:1400]

    display_if_vscode(-devi[300:850, 55:1850])
    display_if_vscode(topo[600:1050, 355:765])
    sleep(1)
    topoc = RGB{N0f8}.(BitmapMaps.__topo_relief(Float32.(z_smooth), CartesianIndices(size(z)), 1))
    display_if_vscode(topoc[600:1050, 355:765])
    sleep(1)
    display_if_vscode(devi[600:1050, 355:765])
    sleep(1)
# => Mostly ok at w = 15

bhc = bh[300:1000, 300:1000]

for window in 3:2:13
    z_smooth = smooth_surface_fir(zc; w = window)
    img = RGB{Float32}.(bhc, # red
            zc,           # green
            Float32.(z_smooth)) # blue
    # Test contours
    conto = __elev_contours(img, minlen, vthick, vdist)
    @show window
    display_if_vscode(conto)
    sleep(1)
end
# => Good at 13



#
# Extract "forest amplitude"
#


# Replace is_forest with bumb_height
forest_cells_min = 20000 #9062
w = 69
cutoff⁻ = 1.5
cutoff⁺ =  4.6477094f0
z_max = 680.0f0
b = bumpiness(z; w, cutoff⁻, z_max, cutoff⁺) 
# Forests are interspersed with high and low values. Smooth that out.
d = smooth_surface_fir(b; w = 15)
display_if_vscode(d[600:850, 500:600])
b[600:1100, 600:1100]
d[600:1100, 600:1100]
# Build a Segmented Image, each forest has its own label.
# Each forest must be large enough to exclude typical rough cliffs (long and slender shapes wihtout much area).
# One (or more) of the segments cover the non-foresty areas.


g = felzenszwalb(d, 2, forest_cells_min)
# Convert segmented image to a normal image matrix
map(i-> Gray{Bool}(segment_mean(g, i) > 0.15), labels_map(g))

for seg in 10:20
 @show segment_mean(g, seg)
end

for k in 2.0:0.5:11.0
    g = felzenszwalb(d, k, forest_cells_min)
    display_if_vscode(map(i -> get_random_color(i), labels_map(g)))
    @show k length(g.segment_labels)
end
# => felzenszwalb gives close to unacceptable results here

for k in 0.75:0.01:0.85
    g = unseeded_region_growing(d, k)
    display_if_vscode(map(i -> get_random_color(i), labels_map(g)))
    @show k length(g.segment_labels)
end
# => unseeded_region_growing easily ends up with strange divisions
#    The result in sens




for k in 0.79:0.1:1.85
    g = region_splitting(d, func_homogeneous(k))
    display_if_vscode(map(i -> get_random_color(i), labels_map(g)))
    @show k length(g.segment_labels)
end
# => region_splitting gives unacceptable results here


for k in 0.15:0.01:0.29
    g = fast_scanning(d, k)
    f_to_be_removed = i -> (segment_pixel_count(g, i) < forest_cells_min)
    f_pick_merge_partner = (i,j) -> (-segment_pixel_count(g, j))
    h = prune_segments(g, to_be_removed, f_pick_merge_partner)
    display_if_vscode(map(i -> get_random_color(i), labels_map(h)))
    @show k length(h.segment_labels)
end
# => fast scanning & pruning gives good results for k = 0.2
#seg = prune_segments(seg, i->(segment_pixel_count(seg,i)<50), (i,j)->(-segment_pixel_count(seg,j)))

g = fast_scanning(d, 0.2)
f_to_be_removed = i -> (segment_pixel_count(g, i) < forest_cells_min)
f_diff(rem_label, neigh_label) = segment_mean(g, neigh_label)

h = prune_segments(g, f_to_be_removed, f_diff)
map(i -> segment_mean(h, i), labels_map(h))
# => Works well, put this in 




RGB{N0f8}.(BitmapMaps.__topo_relief(Float32.(z .- Δz_ideal), CartesianIndices(size(z)), 1))
devi = Δz .- Δz_ideal
@show extrema(Float32.(devi))

for level in 2.0:1:650
    bwimg = map(z -> Gray{Bool}(Float32(z) > level && Float32(z) < (level + 1)), z_smooth)
    display_if_vscode(bwimg[1:1000, 1:1000])
    sleep(0.01)
end
=#