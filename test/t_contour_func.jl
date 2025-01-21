# This optimises parameter and algorithms for identifying forest areas.
#
# This is best done in a slightly large area, which we won't 
# include in the repo resources. Instead, we rely on local files,
# which can be downloaded without creating a user account.
using Test
using BitmapMaps
using BitmapMaps: fir_lp_coefficients, centered, open_as_temp_in_gimp
using BitmapMaps: imfilter, mapwindow, mapwindow!, define_builder, get_random_color
using BitmapMaps: _elev_contours, __elev_contours, Gray, load, felzenszwalb
using BitmapMaps: contour_lines_overlay, CONTOUR_FNAM, CONSOLIDATED_FNAM, TOPORELIEF_FNAM, COMPOSITE_FNAM
using BitmapMaps: RGB, RGBA, labels_map, N0f8, segment_mean, segment_pixel_count, dilate!
using BitmapMaps: CompositeDestinationOver, cartesian_index_string, imfilter, bumpy_patch
import BitmapMaps: bumpy_patch, bumpiness, blue, red, green, smooth_surface_fir, roughness_of_surface_fir
import ImageSegmentation
using ImageSegmentation: unseeded_region_growing, region_splitting, fast_scanning, prune_segments
using Statistics: quantile, mean
using Random


function find_minimimzing_param(f, range)
    val, index = findmin(f_evaluate_fir, range)
    val, range[index]
end

#
# A test image
#
function mountain(r)
    a = 650 / 2
    λ = 1000√2
    χ  =  π * (2r / λ - 1)
    a * (1 - (tanh(χ) ))
end
function mountain(P::Matrix)
    R = CartesianIndices(size(P))
    M = zero.(P)
    for I in R
        xo, yo = I.I
        xc = xo - 1000
        yc = yo - 1000
        r = hypot(xc, yc)
        z = mountain(r)
        M[I] = z
    end
    M
end

function forest(r, θ)
    a = 4
    λ = 5
    χ = r / λ
    va = 0.5 + 0.5 * sin(2 * π * r * θ / λ)
    Float32(a * (1 + va * sin(2 * π * χ)))
end
function forest(P)
    R = CartesianIndices(size(P))
    M = zero.(P)
    for I in R
        xo, yo = I.I
        xc = xo - 1000
        yc = yo - 1000
        r = hypot(xc, yc)
        θ = atan(yc, xc)
        #
        # Add forest in some sectors
        if mod(θ - π / 12, π / 3) > π / 6
            # Wave above 0
            z = forest(r, θ)
        else
            z = 0.0f0
        end
        M[I] = z
    end
    M
end
function noise(P)
    Random.seed!(1) # For consistency
    R = CartesianIndices(size(P))
    M = zero.(P)
    for I in R
        xo, yo = I.I
        xc = xo - 1000
        yc = yo - 1000
        θ = atan(yc, xc)
        # Add noise in some sectors
        if mod(θ, π / 3) > π / 6
            # Wave above 0
            z = Float32(0.5 * randn())
        else
            z = 0.0f0
        end
        M[I] = z
    end
    M
end
function test_terrain(;P = Gray{Float32}[rand()  for i in 1:2000, j in 1:2000])
    mountain(P) .+ forest(P) .+ noise(P)
end

