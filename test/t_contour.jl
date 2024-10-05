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
using BitmapMaps: CompositeDestinationOver, cartesian_index_string, imfilter
import BitmapMaps: is_forest, bumpiness, blue, red, green, smoothed_surface_fir, roughness_of_surface_fir
using Statistics: quantile
using Random



smb = define_builder(;pth = "BitmapMaps/1 1 18145 6894028 71625 6971028", 
    nrc =(3,2), 
    cell_to_utm_factor = 1, 
    southwest_corner = (48705, 6938028),
    density_pt_m⁻¹ = 11811
    )
for i in 1:6
    @test ispath(full_folder_path(smb[i])) # If not:  run_bitmapmap_pipeline(smb)
end

# Separating forest and houses from ridges is pretty hard. We need to use the full 
# resolution (cell_to_utm_factor = 1) to do it well.

function show_results_of_parameters(z, bcg, w)
    b = bumpiness(z; w, cutoff⁺ = 10)
    cou = count(x -> x == 1.0f0, b)
    bmps = filter(x -> x > 0.0f0, b)
    cou = length(bmps)
    mea = Float32(sum(ρ -> ρ / length(bmps), bmps))
    mcou = count(ρ -> ρ == 1.0f0, bmps)
    qua99 = Float32(quantile(bmps, 0.99))
    printstyled(lpad("w = $w ", 15), 
        lpad("count  = $cou ", 20), 
        lpad("maxcount  = $mcou ", 20), 
        lpad("mean(bumps) = $mea ", 30), 
        lpad("99quantile(bumps) = $qua99  ", 30),"\n", color = :green)
    b_green = map(ρ -> RGBA{N0f8}(Float32(ρ == 1.0f0), Float32(ρ^2), Float32(ρ == 1.0f0), Float32(ρ > 0.0f0)), b)
    c = CompositeDestinationOver.(b_green, 0.2 .* bcg)
    display_if_vscode(c)
    c
end

shno = 5
bcg = load(joinpath(full_folder_path(smb[shno]), COMPOSITE_FNAM))
z =  transpose(readclose(joinpath(full_folder_path(smb[shno]), CONSOLIDATED_FNAM)).A[:, :, 1]);
for w in 61:2:81
    show_results_of_parameters(z, bcg, w)
end
# Results minimizing 'count'

# 1 w = 69     count  = 510339       maxcount  = 0       mean(bumps) = 0.2349821 99quantile(bumps) = 0.52533156
# 2 w = 67    count  = 1440404       maxcount  = 0      mean(bumps) = 0.26742142 99quantile(bumps) = 0.6453081 
# 3 w = 67     count  = 358744       maxcount  = 0       mean(bumps) = 0.2559296  99quantile(bumps) = 0.59335  
# 4 w = 73      count  = 10023       maxcount  = 0      mean(bumps) = 0.25086126 99quantile(bumps) = 0.65424484 
# 5 w = 67     count  = 668411       maxcount  = 0       mean(bumps) = 0.2725699 99quantile(bumps) = 0.6557525 
# 6 w = 69     count  = 703201       maxcount  = 0       mean(bumps) = 0.2810644 99quantile(bumps) = 0.6924516  

c = show_results_of_parameters(z, bcg, 67)

save_png_with_phys(joinpath(full_folder_path(smb[shno]), "filterfind.png"), c)

gimp_path = "C:\\Program Files\\GIMP 2\\bin\\gimp-2.10.exe"
for shno in 1:6
    fnam = joinpath(full_folder_path(smb[shno]), "filterfind.png")
    @show fnam
    @assert isfile(fnam)
    cmd = `$gimp_path "$fnam"`
    @show cmd
    @async run(cmd)
end


#
# NEXT ITERATION for finding parameters
#
# It seems w = 65-69 is a good range for the foresty sheets.
# It seems that cap⁺ = 0.6228 * 10 = 6.228f0  should cover 99% of actual forest. 
# Let's revise this function, and also cut values above max instead of capping.

# Now that we don't primarily want to emphasize the cut-off pixels, we change 
# the exponent so as to show the extent of areas better.
# Let's revise this function.
function show_results_of_parameters(z, bcg, w)
    b = bumpiness(z; w)
    cou = count(x -> x == 1.0f0, b)
    bmps = filter(x -> x > 0.0f0, b)
    cou = length(bmps)
    mea = Float32(sum(ρ -> ρ / length(bmps), bmps))
    mcou = count(ρ -> ρ == 1.0f0, bmps)
    qua99 = Float32(quantile(bmps, 0.99))
    printstyled(lpad("w = $w ", 15), 
        lpad("count  = $cou ", 20), 
        lpad("maxcount  = $mcou ", 20), 
        lpad("mean(bumps) = $mea ", 30), 
        lpad("99quantile(bumps) = $qua99  ", 30),"\n", color = :green)
    b_green = map(ρ -> RGBA{N0f8}(Float32(ρ == 1.0f0), Float32(ρ^0.2), Float32(ρ == 1.0f0), Float32(ρ > 0.0f0)), b)
    c = CompositeDestinationOver.(b_green, 0.2 .* bcg)
    display_if_vscode(c)
    c
end




shno = 1
bcg = load(joinpath(full_folder_path(smb[shno]), COMPOSITE_FNAM))
z =  transpose(readclose(joinpath(full_folder_path(smb[shno]), CONSOLIDATED_FNAM)).A[:, :, 1]);
for w in 65:2:73
    show_results_of_parameters(z, bcg, w)
end

# Results minimizing 'count'
# 1         w = 69     count  = 165360       maxcount  = 0      mean(bumps) = 0.52542007 99quantile(bumps) = 0.90849406 
# 2         w = 67     count  = 623331       maxcount  = 0      mean(bumps) = 0.55886495 99quantile(bumps) = 0.9517818  
# 3         w = 65     count  = 143469       maxcount  = 0       mean(bumps) = 0.5439831 99quantile(bumps) = 0.9308489 
# 4         w = 71     count  = 116181       maxcount  = 0      mean(bumps) = 0.57334054 99quantile(bumps) = 0.97327125
# 5         w = 69     count  = 305137       maxcount  = 0      mean(bumps) = 0.56516373 99quantile(bumps) = 0.9562366  
# 6         w = 67     count  = 332824       maxcount  = 0       mean(bumps) = 0.5703782 99quantile(bumps) = 0.9618552 

c = show_results_of_parameters(z, bcg, 67)
fnam = joinpath(full_folder_path(smb[shno]), "filterfind.png")
save_png_with_phys(fnam, c)
@assert isfile(fnam)
@async run(`$gimp_path "$fnam"`)



#
# NEXT ITERATION for finding parameters
#
# It seems w = 69 is a good value for the foresty sheets.
# It seems that cutoff⁺ = 0.90849406 * 6.228f0 = 5.658101f0  should cover 99% of actual forest. 
# To really focus on actual forest, we lower z_max to a realistic local forest elevation boundary.
# Let's revise this function.
# Let's revise this function.
function show_results_of_parameters(z, bcg, cutoff⁺)
    b = bumpiness(z; cutoff⁺)
    cou = count(x -> x == 1.0f0, b)
    bmps = filter(x -> x > 0.0f0, b)
    cou = length(bmps)
    mea = Float32(sum(ρ -> ρ / length(bmps), bmps))
    mcou = count(ρ -> ρ == 1.0f0, bmps)
    qua99 = Float32(quantile(bmps, 0.99))
    printstyled(lpad("count  = $cou ", 20), 
        lpad("mean(bumps) = $mea ", 30), 
        lpad("99quantile(bumps) = $qua99  ", 30),"\n", color = :green)
    b_green = map(ρ -> RGBA{N0f8}(Float32(ρ == 1.0f0), Float32(ρ^0.2), Float32(ρ == 1.0f0), Float32(ρ > 0.0f0)), b)
    c = CompositeDestinationOver.(b_green, 0.2 .* bcg)
    display_if_vscode(c)
    c
end


shno = 1
bcg = load(joinpath(full_folder_path(smb[shno]), COMPOSITE_FNAM))
z =  transpose(readclose(joinpath(full_folder_path(smb[shno]), CONSOLIDATED_FNAM)).A[:, :, 1]);
show_results_of_parameters(z, bcg, 5.658101f0) 
# count  = 155611      mean(bumps) = 0.57326555 99quantile(bumps) = 0.93999803
show_results_of_parameters(z, bcg, 5.658101f0 * 0.93999803) 
# count  = 154054       mean(bumps) = 0.6056124 99quantile(bumps) = 0.9569906
show_results_of_parameters(z, bcg, 5.658101f0 * 0.93999803 * 0.957) 
# count  = 152514       mean(bumps) = 0.6289043 99quantile(bumps) = 0.9642101 
show_results_of_parameters(z, bcg, 5.658101f0 * 0.93999803 * 0.957 * 0.9642)
# count  = 150988        mean(bumps) = 0.648564 99quantile(bumps) = 0.9715482  
show_results_of_parameters(z, bcg, 5.658101f0 * 0.93999803 * 0.957 * 0.9642 * 0.9715482) 
# count  = 149478      mean(bumps) = 0.66405535 99quantile(bumps) = 0.9747605 
# Let's decide on cutoff⁺ =  5.658101f0 * 0.93999803 * 0.957 * 0.9642 * 0.9715482 * 0.9747605  = 4.6477094f0 Revising:

# Let's revise this function.
# New parameter for optimization, and use linear bumps.
function show_results_of_parameters(z, bcg, cutoff⁻)
    b = bumpiness(z; cutoff⁻)
    cou = count(x -> x == 1.0f0, b)
    bmps = filter(x -> x > 0.0f0, b)
    cou = count(x -> x == 1.0f0, b)
    bmps = filter(x -> x > 0.0f0, b)
    cou = length(bmps)
    mea = Float32(sum(ρ -> ρ / length(bmps), bmps))
    mcou = count(ρ -> ρ == 1.0f0, bmps)
    qua99 = Float32(quantile(bmps, 0.99))
    printstyled("# ", lpad("count  = $cou ", 20), 
        lpad("mean(bumps) = $mea ", 30), 
        lpad("99quantile(bumps) = $qua99  ", 30),"\n", color = :green)
    b_green = map(ρ -> RGBA{N0f8}(Float32(ρ == 1.0f0), Float32(ρ), Float32(ρ == 1.0f0), Float32(ρ > 0.0f0)), b)
    c = CompositeDestinationOver.(b_green, 0.4 .* bcg)
    display_if_vscode(c)
    c
end


show_results_of_parameters(z, bcg, 2.5)
#     count  = 149353       mean(bumps) = 0.6778929 99quantile(bumps) = 0.97844714  
show_results_of_parameters(z, bcg, 2.0)
#     count  = 278216       mean(bumps) = 0.5858495 99quantile(bumps) = 0.9621004 
show_results_of_parameters(z, bcg, 1.5)
#     count  = 499945      mean(bumps) = 0.49092907 99quantile(bumps) = 0.93621135 
show_results_of_parameters(z, bcg, 1.0)
#     count  = 874794      mean(bumps) = 0.39383477 99quantile(bumps) = 0.89969987
shno = 5
bcg = load(joinpath(full_folder_path(smb[shno]), COMPOSITE_FNAM))
z =  transpose(readclose(joinpath(full_folder_path(smb[shno]), CONSOLIDATED_FNAM)).A[:, :, 1]);
c = show_results_of_parameters(z, bcg, 1.5)
c[CartesianIndices((100:800, 1:500))]
c[CartesianIndices((1000:1800, 500:1200))]
# Let's change the lower cutoff default to 1.5 (and the tree-line to 680)

# Check these parameters results in Gimp:
for shno = 1:6
    @show shno
    bcg = load(joinpath(full_folder_path(smb[shno]), COMPOSITE_FNAM))
    z =  transpose(readclose(joinpath(full_folder_path(smb[shno]), CONSOLIDATED_FNAM)).A[:, :, 1]);
    b = bumpiness(z)
    b_green = map(ρ -> RGBA{N0f8}(0, Float32(ρ), 0, Float32(ρ > 0.0f0)), b)
    c = CompositeDestinationOver.(b_green, 0.4 .* bcg)
    fnam = joinpath(full_folder_path(smb[shno]), "filterfind.png")
    save_png_with_phys(fnam, c)
    @assert isfile(fnam)
    @async run(`$gimp_path "$fnam"`)
end
 



# Let's check with a more built-up area:
smb = define_builder(;pth = "BitmapMaps/1 1 18145 6894028 71625 6971028", 
    nrc =(3,3), 
    cell_to_utm_factor = 1, 
    southwest_corner = (34755, 6921336),
    density_pt_m⁻¹ = 11811
    )

for i in 1:9
    @test ispath(full_folder_path(smb[i])) # If not:  run_bitmapmap_pipeline(smb)
end

# Check these parameters results in Gimp (which should be opened by user, or else
# wont' be visible):
for shno in 1:9
    println(cartesian_index_string(smb, shno), " of ", cartesian_index_string(smb))
    bcg = load(joinpath(full_folder_path(smb[shno]), COMPOSITE_FNAM))
    z =  transpose(readclose(joinpath(full_folder_path(smb[shno]), CONSOLIDATED_FNAM)).A[:, :, 1]);
    b = bumpiness(z)
    b_green = map(ρ -> RGBA{N0f8}(0, Float32(ρ), 0, Float32(ρ > 0.0f0)), b)
    c = CompositeDestinationOver.(b_green, 0.4 .* bcg)
    fnam = joinpath(full_folder_path(smb[shno]), "filterfind.png")
    save_png_with_phys(fnam, c)
    @assert isfile(fnam)
    cmd = `$gimp_path "$fnam"`
    @show cmd
    @async run(cmd)
end





# This looks ready for some Gaussian blurring (which is also a low-pass filter).
shno = 5
println(cartesian_index_string(smb, shno), " of ", cartesian_index_string(smb))
bcg = load(joinpath(full_folder_path(smb[shno]), COMPOSITE_FNAM))
z =  transpose(readclose(joinpath(full_folder_path(smb[shno]), CONSOLIDATED_FNAM)).A[:, :, 1]);
b = bumpiness(z)
b[CartesianIndices((1000:1600, 500:1200))]
c = map(β -> Gray{N0f8}(β > 0.0), b)
c[CartesianIndices((1000:1600, 500:1200))]
d = imfilter(b, Kernel.gaussian(2))
d[CartesianIndices((1000:1600, 500:1200))]
e = map(β -> Gray{N0f8}(β > 0.0), d)
e[CartesianIndices((1000:1600, 500:1200))]
f = CompositeDestinationOver.(map(ρ -> RGBA{N0f8}(0, Float32(ρ), 0, Float32(0.2 * (ρ > 0))), e), 0.9 .* bcg)
f[CartesianIndices((2500:3000, 1000:1800))]
# Forest is applied on a few cliffs below 680 m. Let's do away with the smaller segments.
forest_pixels_min = 9062#1750 # Iterated to minimum
g = felzenszwalb(e, 1, forest_pixels_min)
h = map(i->segment_mean(g, i), labels_map(g))
i = CompositeDestinationOver.(map(ρ -> RGBA{N0f8}(0, Float32(ρ), 0, Float32(0.3 * (ρ > 0))), h), 1.0 .* bcg)

i[CartesianIndices((2500:3000, 1:800))]

#
# A test image
#


function mountain(r)
    a = 650 / 2
    λ = 1000√2
    χ  =  π * (2r / λ - 1)
    a * (1 - (tanh(χ) ))
end
function forest(r, θ)
    a = 4
    λ = 5
    χ = r / λ
    va = 0.5 + 0.5 * sin(2 * π * r * θ / λ)
    a * (1 + va * sin(2 * π * χ))
end


function test_elevations()
    P = Gray{Float32}[rand()  for i in 1:2000, j in 1:2000]
    R = CartesianIndices(size(P))
    for I in R
        xo, yo = I.I
        xc = xo - 1000
        yc = yo - 1000
        r = hypot(xc, yc)
        θ = atan(yc, xc)
        #
        z = mountain(r)
        # Add forest in some sectors
        if mod(θ - π / 12, π / 3) > π / 6
            # Wave above 0
            z += forest(r, θ)
        end
        # Add noise in some
        if mod(θ, π / 3) > π / 6
            # Wave above 0
            z += 0.5 *randn()
        end
        P[I] = z
    end
    P
end


z = test_elevations()
# topo relief
topo = RGB{N0f8}.(BitmapMaps.__topo_relief(Float32.(z), CartesianIndices(size(z)), 1))
display_if_vscode(topo[900:1050, 755:765])

vthick = [1, 3, 5]
vdist = [20, 100, 1000] # Contour spacings.
minlen = 6
# Pre-process as in _elev_contours (won't be visible)
isfo = Float32.(is_forest(z))
display_if_vscode(isfo)
z_smooth = smoothed_surface_fir(z; w = 13)
display_if_vscode(z_smooth)

for level in 2.0:1:650
    bwimg = map(z -> Gray{Bool}(Float32(z) > level && Float32(z) < (level + 1)), z_smooth)
    display_if_vscode(bwimg[1:1000, 1:1000])
    sleep(0.01)
end


img = RGB{Float32}.(isfo, # red
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


# DELETE

# The low-pass filter deforms the large-scale terrain a lot. We must do better.

zc = z[300:1000, 300:1000]
RC = CartesianIndices(size(zc))
for window in 3:2:27
    z_smooth = smoothed_surface_fir(zc; w = window)
    topo = RGB{N0f8}.(BitmapMaps.__topo_relief(Float32.(z_smooth), RC, 1))
    @show window
    display_if_vscode(topo)
end
# => Mostly ok at 15

isfoc = isfo[300:1000, 300:1000]

for window in 3:2:13
    z_smooth = smoothed_surface_fir(zc; w = window)
    img = RGB{Float32}.(isfoc, # red
            zc,           # green
            Float32.(z_smooth)) # blue
    # Test contours
    conto = __elev_contours(img, minlen, vthick, vdist)
    @show window
    display_if_vscode(conto)
    sleep(1)
end
# => Good at 13