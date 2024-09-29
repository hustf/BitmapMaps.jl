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
import BitmapMaps: Kernel
using BitmapMaps.Kernel: gaussian
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
# resolution (cell_to_utm_factor = 1) to have much chance.

function bumpiness(z; w = 67, cutoff⁻ = 2.5, z_max = 6000.0f0, cap⁺ = 16.0f0)
    isodd(w) || throw(ArgumentError("w must be odd, not $w"))
    # Prepare coefficients for wide high-pass filter. 
    #     Coefficients, including a window.
    c = Float32.(BitmapMaps.fir_hp_coefficients(w))
    # The local bumps (+ and -)
    r = imfilter(z, (c, transpose(c)), BitmapMaps.FIRTiled())
    r1 = r .* (z .< z_max)
    # Take absolute value, and cap the large values, which may be artifacts, 
    # and would not affect significant areas anyway.
    # NOTE: Maybe drop the large values altogether?
    map(r1) do ρ
        mag = abs(ρ)
        mag < cutoff⁻ ? Gray{Float32}(0.0f0) : mag > cap⁺ ? Gray{Float32}(1.0f0) : Gray{Float32}(mag / cap⁺)
    end
end



function show_results_of_parameters(z, bcg, w)
    b = bumpiness(z; w, cap⁺ = 10)
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

shno = 6
bcg = load(joinpath(full_folder_path(smb[shno]), COMPOSITE_FNAM))
z =  transpose(readclose(joinpath(full_folder_path(smb[shno]), CONSOLIDATED_FNAM)).A[:, :, 1]);
for w in 61:2:81
    show_results_of_parameters(z, bcg, w)
end
# Results minimizing 'count'
# shno = 1      w = 69     count  = 167032      maxcount  = 21       mean(bumps) = 0.3310675 99quantile(bumps) = 0.6228346 
# shno = 2      w = 65     count  = 641551     maxcount  = 331      mean(bumps) = 0.35834447 99quantile(bumps) = 0.7264333
# shno = 3      w = 65     count  = 146164      maxcount  = 99      mean(bumps) = 0.34591758 99quantile(bumps) = 0.68322724 
# shno = 4      w = 71     count  = 128875    maxcount  = 2652      mean(bumps) = 0.40135497      99quantile(bumps) = 1.0 
# shno = 5      w = 69     count  = 315338     maxcount  = 569       mean(bumps) = 0.3643127 99quantile(bumps) = 0.7528747 
#        6       w = 67     count  = 346976     maxcount  = 728      mean(bumps) = 0.37078118 99quantile(bumps) = 0.78817546 

c = show_results_of_parameters(z, bcg, 67)
c[CartesianIndices((2100:size(z, 1), 1000:size(z, 2)))]
c[CartesianIndices((2100:2500, 1000:1400))]

save_png_with_phys(joinpath(full_folder_path(smb[shno]), "filterfind.png"), c)

# Gimp should be opened by user first, or else
# wont' be visible:
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
function bumpiness(z; w = 67, cutoff⁻ = 2.5, z_max = 6000.0f0, cutoff⁺ = 6.228f0)
    isodd(w) || throw(ArgumentError("w must be odd, not $w"))
    # Prepare coefficients for wide high-pass filter. 
    #     Coefficients, including a window.
    c = Float32.(BitmapMaps.fir_hp_coefficients(w))
    # The local bumps (+ and -)
    r = imfilter(z, (c, transpose(c)), BitmapMaps.FIRTiled())
    r1 = r .* (z .< z_max)
    # Take absolute value, and drop the large values, which are most likely artifacts
    # or houses. 
    map(r1) do ρ
        mag = abs(ρ)
        mag < cutoff⁻ ? Gray{Float32}(0.0f0) : mag > cutoff⁺ ? Gray{Float32}(0.0f0) : Gray{Float32}(mag / cutoff⁺)
    end
end

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
function bumpiness(z; w = 69, cutoff⁻ = 2.5, z_max = 600.0f0, cutoff⁺ = 5.658101f0)
    isodd(w) || throw(ArgumentError("w must be odd, not $w"))
    # Prepare coefficients for wide high-pass filter. 
    #     Coefficients, including a window.
    c = Float32.(BitmapMaps.fir_hp_coefficients(w))
    # The local bumps (+ and -)
    r = imfilter(z, (c, transpose(c)), BitmapMaps.FIRTiled())
    r1 = r .* (z .< z_max)
    # Take absolute value, and drop the large values, which are most likely artifacts
    # or houses. 
    map(r1) do ρ
        mag = abs(ρ)
        mag < cutoff⁻ ? Gray{Float32}(0.0f0) : mag > cutoff⁺ ? Gray{Float32}(0.0f0) : Gray{Float32}(mag / cutoff⁺)
    end
end
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
function bumpiness(z; w = 69, cutoff⁻ = 2.5, z_max = 700.0f0, cutoff⁺ = 4.6477094f0)
    isodd(w) || throw(ArgumentError("w must be odd, not $w"))
    # Prepare coefficients for wide high-pass filter. 
    #     Coefficients, including a window.
    c = Float32.(BitmapMaps.fir_hp_coefficients(w))
    # The local bumps (+ and -)
    r = imfilter(z, (c, transpose(c)), BitmapMaps.FIRTiled())
    r1 = r .* (z .< z_max)
    # Take absolute value, and drop the large values, which are most likely artifacts
    # or houses. 
    map(r1) do ρ
        mag = abs(ρ)
        mag < cutoff⁻ ? Gray{Float32}(0.0f0) : mag > cutoff⁺ ? Gray{Float32}(0.0f0) : Gray{Float32}(mag / cutoff⁺)
    end
end

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
function bumpiness(z; w = 69, cutoff⁻ = 1.5, z_max = 680.0f0, cutoff⁺ =  4.6477094f0)
    isodd(w) || throw(ArgumentError("w must be odd, not $w"))
    # Prepare coefficients for wide high-pass filter. 
    #     Coefficients, including a window.
    c = Float32.(BitmapMaps.fir_hp_coefficients(w))
    # The local bumps (+ and -)
    r = imfilter(z, (c, transpose(c)), BitmapMaps.FIRTiled())
    r1 = r .* (z .< z_max)
    # Take absolute value, and drop the large values, which are most likely artifacts
    # or houses. 
    map(r1) do ρ
        mag = abs(ρ)
        mag < cutoff⁻ ? Gray{Float32}(0.0f0) : mag > cutoff⁺ ? Gray{Float32}(0.0f0) : Gray{Float32}(mag / cutoff⁺)
    end
end

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

"""
    is_forest(z; forest_cell_min = 9062, w = 69, cutoff⁻ = 1.5, cutoff⁺ =  4.6477094f0, z_max = 680.0f0)
    ---> Matrix{Gray{Bool}}

Used by 'ridges' and 'contours'.

# Arguments

z                 Elevations. The other default parameters assumes a grid spacing of 1 utm meter.
forest_cell_min   Remove smaller forests
w                 Window length for high-pass filter
cutoff⁻           Drop filtered values below
cutoff⁺           Drop filtered values above (often steep ridges or artifacts)
z_max             Drop forest above this elevation
"""
function is_forest(z; forest_cell_min = 9062, w = 69, cutoff⁻ = 1.5, cutoff⁺ =  4.6477094f0, z_max = 680.0f0)
    b = bumpiness(z; w, cutoff⁻, z_max, cutoff⁺)
    d = imfilter(b, Kernel.gaussian(2))
    e = map(β -> Gray{N0f8}(β > 0.0), d)
    g = felzenszwalb(e, 1, forest_pixels_min)
    map(i-> Gray{Bool}(round(segment_mean(g, i))), labels_map(g))
end