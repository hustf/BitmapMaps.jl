# Internals to _elev_contour.
# This is best checked in a slightly large area, which we won't 
# include in the repo resources. Instead, we rely on local files,
# which can be downloaded without creating a user account.
using Test
using BitmapMaps
using BitmapMaps: fir_lp_coefficients, centered
using BitmapMaps: imfilter, mapwindow, mapwindow!
using ImageCore: N0f32, N0f16, N0f8, Normed, scaleminmax, RGB, channelview, GrayA, RGBA, gray, alpha
using ImageCore: red, green, blue, Gray, mosaic
using ImageFiltering: FIRTiled, imfilter!, imgradients, KernelFactors, Fill
import ColorBlendModes
using ColorBlendModes: BlendLighten
import ImageMorphology: erode, dilate


fofo = joinpath(homedir(), "bitmapmaps\\render\\1 1  47675 6929520  54041 6938686")
g = readclose(joinpath(fofo, BitmapMaps.CONSOLIDATED_FNAM))
source_indices = (1:2:size(g)[1], 1:2:size(g)[2])
csi = CartesianIndices(source_indices)
# Make an array of elevation (red) plus filtered steepness (blue).
# The blue channel is added below.
zs = GrayA{Float32}.(g.A[:, :, 1], Array{Float32}(undef, size(g.A[:, :, 1])...))
# This emphasizes small value differences. Inspect filtering with it.
# view(zsc, 7:1400, 7:1400) |> img -> scaleminmax(RGB, extrema(red.(img))...).(img)

# Vector filter coefficients, for finding a good steepness value in e.g. forests.
w = 49 # Window size for low-pass-filter
c = Float32.(fir_lp_coefficients(w))
r = channelview(zs)[1, :, :]
zf = imfilter(r, (c, transpose(c)), FIRTiled());
zfx, zfy  = imgradients(zf, KernelFactors.ando5, "replicate")
steepness = Float32.(hypot.(zfx, zfy))
transpose(reinterpret(Gray{Float32}, steepness))
channelview(zs)[2, :, :] .= steepness
α = channelview(zs)[2, :, :]
side = 1000
no = 200
insp = CartesianIndices((1:side, (9166 - no - side):(9166 - no)))
mimag = extrema(view(α, insp))
mimar = extrema(view(r, insp))
# zoom in
mosaic(view(zf, insp) .|> scaleminmax(mimar...) .|> Gray |> transpose, view(α, insp) .|> scaleminmax(mimag...) .|> Gray |> transpose)
# all
mosaic(zf .|> scaleminmax(extrema(r)...) .|> Gray |> transpose, α .|> scaleminmax(extrema(α)...) .|> Gray |> transpose)



fc = BitmapMaps.func_elev_contour(20f0)
# pre-allocate
bw = Array{Gray{Bool}}(undef, size(csi)...)
# 0.233s
@time mapwindow!(fc, bw, zs, (1, 1); indices = source_indices)

bw = BitmapMaps.strokeify(bw, 5)



imcol = map(bw) do pix
    pix == true && return RGBA{N0f8}(0.173, 0.192, 0.255, 1)
    RGBA{N0f8}(0, 0, 0, 0)
end

ffna = joinpath(fofo, BitmapMaps.CONTOUR_FNAM)
save_png_with_phys(ffna, transpose(imcol))

#
# Combine black-white images:
#
img =     Gray{Bool}[0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0]
overlay = Gray{Bool}[0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 1]
fa!(img, overlay) = reinterpret(Bool, img) .|= reinterpret(Bool, overlay)
# 30.522 ns (0 allocations: 0 bytes)
#@btime fa!($img, $overlay)
fa!(img, overlay)
@test img == Gray{Bool}[0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]

# Couldn't do this with map!
function foi(img, overlay)
    map(zip(img, overlay)) do (x, y)
        if x == true || y == true 
            Gray{Bool}(true)
        else
            Gray{Bool}(false)
        end
    end
end
#  39.516 ns (1 allocation: 80 bytes)
#@btime img = foi($img, $overlay)
img = foi(img, overlay)
@test img == Gray{Bool}[0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]

img == Gray{Bool}[0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
# Canonical way is best way:
fii!(img, overlay) = map!(BlendLighten, img, img, overlay)
# 4.200 ns (0 allocations: 0 bytes)
#@btime fii!($img, $overlay)



##########################
# Removing isolated pixels
##########################

@test overlay == Gray{Bool}[0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 1]

overlay = Gray{Bool}[0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 1]
BitmapMaps.remove_isolated_pixels!(overlay)
@test sum(overlay) == 0

overlay = Gray{Bool}[0 0 0 0
                     0 1 0 0
                     0 0 1 0
                     0 0 0 1]
BitmapMaps.remove_isolated_pixels!(overlay)
@test sum(overlay) == 3 

##########
# Removing islands of pixels
##########
using BitmapMaps: felzenszwalb, segment_pixel_count, labels_map, segment_mean
overlay = Gray{Bool}[0 0 0 0
                     0 1 0 0
                     0 0 1 0
                     0 0 0 1]

img = map(overlay) do pix
    pix == true && return RGB{N0f8}(1.0, 1.0, 1.0)
    RGB{N0f8}(0., 0, 0)
end


img = overlay

segments = felzenszwalb(reinterpret(Bool, img), 1.0, 2)
is_small(i) = segment_pixel_count(segments, i) <= 3
is_blank(i) = segment_mean(segments, i) < 1
issmall = map(labels_map(segments)) do i # i is a label, representing a set of pixels.
    Gray{Bool}(is_small(i) && is_blank(i))
end
function remove_small_islands(img)
    segments = felzenszwalb(reinterpret(Bool, img), 1.0, 2)
    is_large(i) = segment_pixel_count(segments, i) > 3
    is_nonblank(i) = segment_mean(segments, i) == 1
    map(labels_map(segments)) do i # i is a label in segments, representing a set of pixels.
        Gray{Bool}(is_large(i) && is_nonblank(i))
    end
end

img = Gray{Bool}[0 0 0 0
                 0 1 0 0
                 0 1 1 0
                 0 0 1 1]
@test remove_small_islands(img) == img
img = Gray{Bool}[0 0 0 0
                 0 1 0 0
                 0 1 1 0
                 0 1 0 0]
@test remove_small_islands(img) == img
img = Gray{Bool}[0 0 0 0
                 0 1 0 0
                 0 1 0 0
                 1 0 0 0]
@test remove_small_islands(img) == zeros(4,4)
img = Gray{Bool}[0 0 0 0
                 0 1 0 0
                 0 0 0 0
                 0 0 0 0]
@test remove_small_islands(img) == zeros(4,4)
img = Gray{Bool}[0 0 0 0
                 0 1 0 0
                 0 0 1 0
                 0 0 0 1]
@test remove_small_islands(img) == zeros(4,4)
img = Gray{Bool}[0 0 1 0
                 0 1 0 0
                 0 0 1 0
                 0 0 0 1]
@test remove_small_islands(img) == img


function remove_small_islands!(img::Array{Gray{Bool},2})
    segments = felzenszwalb(reinterpret(Bool, img), 1.0, 2)
    labels = labels_map(segments)
    # Count the number of pixels in each segment
    segment_counts = segment_pixel_count(segments)
    # Create a boolean array for valid segments (larger than 3 and nonblank)
    valid_segments = Dict(i => (segment_counts[i] > 3 && segment_mean(segments, i) == 1) for i in keys(segment_counts))
    # Update the img array in-place
    for i in eachindex(img)
        img[i] = Gray{Bool}(valid_segments[labels[i]])
    end
    img
end
img = Gray{Bool}[0 0 1 0
                 0 1 0 0
                 0 0 1 0
                 0 0 0 1]
remove_small_islands!(img)
@test img == Gray{Bool}[0 0 1 0
    0 1 0 0
    0 0 1 0
    0 0 0 1]
img = Gray{Bool}[0 0 1 0
                 0 1 0 0
                 0 0 1 0
                 0 0 0 0]
remove_small_islands!(img)
@test img == zeros(4,4)

img = Gray{Bool}[0 0 1 0
                0 1 0 0
                0 0 1 0
                0 0 0 0]
# 2.767 μs (30 allocations: 6.07 KiB)
#@btime remove_small_islands($img)
# 2.344 μs (32 allocations: 5.34 KiB)
#@btime remove_small_islands!($img)


