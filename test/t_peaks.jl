# NOTE:
# This is work in progress. Also related: t_prominence and mark_utils.jl. We are trying to identify and mark 
# peaks from their prominence ("primary factor"), and while to reusing some functions from 'contour.jl'.
# Marking may be delegated to the svg.
using Test
using BitmapMaps
import ColorTypes
using ColorTypes: RGB, N0f8, Gray, GrayA
using ImageCore: red, green, blue, alpha, gray, mosaic, scaleminmax, channelview
using ImageCore: zeroarray, colorview
using BitmapMaps: indices, gaussian_curvature, gaussian_curvature_smoothed
import ImageMorphology
using ImageMorphology: MaxTree, area_opening, area_closing, diameter_opening, diameter_closing
using ImageMorphology: areas, diameters, local_maxima!, local_minima!, extreme_filter, extreme_filter!
using ImageMorphology: StructuringElements.strel_box, label_components
import ImageSegmentation
using ImageSegmentation: labels_map, segment_labels, seeded_region_growing
using ImageView # This often triggers slow responding
ImageView.closeall()  # This stops the slow response state
#using ImageFiltering: mapwindow!, centered

# Extract an interesting area for identifying mountain peaks.
# This is best done in a slightly large area, which we don't
# include in the repo resources. Instead, we rely on local files.

fofo = joinpath(homedir(), "bitmapmaps\\render\\1 1  47675 6929520  54041 6938686")
g = readclose(joinpath(fofo, BitmapMaps.CONSOLIDATED_FNAM))
# Likely, this is the most peaky part
ne_utm = (49790.9, 6933344.6)
sw_utm = (49198.0, 6932656.0)
sw_ind = indices(g, sw_utm)
ne_ind = indices(g, ne_utm)
nw_ind = CartesianIndex(sw_ind[1], ne_ind[2])
se_ind = CartesianIndex(ne_ind[1], sw_ind[2])
cis = nw_ind:se_ind
# Make an RGB array of elevation (red). Other channels will be defined below.
zs = RGB{Float32}.(g[cis], Array{Float32}(undef, size(cis)...), Array{Float32}(undef, size(cis)...))
re = channelview(zs)[1, :, :]
gr = channelview(zs)[2, :, :]
bl = channelview(zs)[3, :, :]

# Inspection
transpose(colorview(RGB, scaleminmax(extrema(re)...).(re), zeroarray, zeroarray))
# Inspect function with color scaling. The window may appear behind others.
ImageView.closeall()
cnv = imshow(zs)["gui"]["canvas"]
sleep(0.5)
function insp(a, b, c)
    img = transpose(colorview(RGB, Float32(a) .* (scaleminmax(extrema(re)...).(re)),
                             Float32(b) .* (scaleminmax(extrema(gr)...).(gr)),
                             Float32(c) .* (scaleminmax(extrema(bl)...).(bl))))
    imshow(cnv, img)
    nothing
end
#
# After preparation: Identify peaks.
#

# This flattens the lower 10% to the maximum value in that area.
#min_area = *(size(re)...) * 0.1
#gr .= area_closing(re, MaxTree(re; rev = true); min_area)
#insp(0,1,0)
# This flattens the upper 10% to the minimum value in that area.
#bl .= area_opening(re, MaxTree(re); min_area) # How does it work? Fast?
#insp(0,0,1)
# Let's extract more info from tree.....
# This contains zero except at local maxima.
# Maxima elements contains a unique id
local_minima!(gr, re, MaxTree(re, rev = true, connectivity = 2))
count(>(0), gr) # 157
insp(0.5,100,0.0) # These are not true minima in the full neighbourhood
local_minima!(gr, re, connectivity = 2)
insp(0.8,1,0)
@test count(>(0), gr) == 157

# We make our own feature identifier and can veryfiy that the first one was correct.
# Inspection through ImageViews.jl (on Windows) is incorrect, because the shown coordinate is not
# where your mouse hovers.
extreme_filter!(min, gr, re, strel_box((3,3)))
gr .= gr .== re
@test count(>(0), gr) == 157
#

bl .= gaussian_curvature(re) .< -0.2
@test sum(bl) == 3
insp(1, 0, 1)
# It shall be easier to inspect / verify the criterion with some contours. We put contours on the green channel.
zs = GrayA{Float32}.(re, Array{Float32}(undef, size(re)...))
channelview(zs)[2, :, :] .= BitmapMaps.terrain_steepness(zs)
BitmapMaps.mapwindow!(BitmapMaps.func_elev_contour(10), gr, zs, (1, 1))
insp(1, 0.1, 1)


# Results were not good. Let's smooth the surface prior to finding saddle points.

bl .= BitmapMaps.gaussian_curvature_smoothed(zf) .< -0.005
@test count(>(0), bl) == 554
insp(1,0,1)
# These values look good. The number of saddle points seem large, but most are grouped. We can reduce the number of elevations to check easily.

function sorted_saddle_elevations(z)
    labelled_saddle_points = label_components(gaussian_curvature_smoothed(z) .< -0.005);
    num_saddle_segments = maximum(labelled_saddle_points)
    # Initialize arrays to store min and max values
    min_values = zeros(Float32, num_saddle_segments)
    max_values = zeros(Float32, num_saddle_segments)
    # Loop through each segment to find min and max values
    for i in 1:num_segments
        segment_values = z[labelled_saddle_points .== i]
        min_values[i] = minimum(segment_values)
        max_values[i] = maximum(segment_values)
    end
    _sorted_saddle_elevations(min_values, max_values)
end

function _sorted_saddle_elevations(min_values, max_values)
    merged_intervals = merge_intervals(min_values, max_values)
    # Extract unique values from the merged intervals
    values = eltype(merged_intervals[1])[]
    for interval in merged_intervals
        append!(values, interval)
    end
    sort(unique(values))
end

function merge_intervals(min_values, max_values; max_gap_between_intervals = 1.0f0)
    intervals = sort(collect(zip(min_values, max_values)))
    merged_intervals = []
    current_interval = intervals[1]
    for i in 2:length(intervals)
        next_interval = intervals[i]
        if next_interval[1] - max_gap_between_intervals <= current_interval[2]
            println("Merged:     $current_interval  $next_interval")
            current_interval = (current_interval[1], max(current_interval[2], next_interval[2]))
        else
            println("Not merged: $current_interval  $next_interval")
            push!(merged_intervals, current_interval)
            current_interval = next_interval
        end
    end
    push!(merged_intervals, current_interval)
    merged_intervals
end

sorted_saddle_elevations(re)


sorted_saddle_elevations(min_values, max_values)



channelview(zs)[2, :, :] .= gaussian_curvature(zf) .< -0.005
channelview(zs)[1, :, :] .= re
saddle_point_elevations = map(zs) do pix
    alpha(pix) == true ? gray(pix) : 0.0f0
end


M = saddle_point_elevations[305:322, 310:338]
M[2,2:4] .= 900.0 .+ rand()

M1 = M[:, 1:12]
println(M1)
display(transpose(Gray.(saddle_point_elevations[305:322, 310:338])))


# Convert the input matrix to an image
input_image = colorview(Gray, saddle_point_elevations)

# Perform image segmentation to label the segments
segmented_image = label_components(input_image .> 0.0)

# Extract the number of segments
num_segments = maximum(segmented_image)

# Initialize arrays to store min and max values
min_values = zeros(Float32, num_segments)
max_values = zeros(Float32, num_segments)

# Loop through each segment to find min and max values
for i in 1:num_segments
    segment_values = re[segmented_image .== i]
    min_values[i] = minimum(segment_values)
    max_values[i] = maximum(segment_values)
end

# Print the results
for i in 1:num_segments
    println("Segment $i: min = ", min_values[i], ", max = ", max_values[i])
end













saddle_point_labels = label_components(gaussian_curvature(zf) .< -0.005)

seeded_region_growing(re, saddle_point_elevations)


labels_map(saddle_point_labels)
segment_labels(saddle_point_labels)

is_saddle_point(I) = gaussian_curvature(zf)[I] .< -0.005

bl .=














mapwindow!(is_local_minimum, gr, re, (3, 3); border = "replicate")

mapwindow!(is_local_minimum, gr, re, ImageMorphology.StructuringElements.strel_box((3,3)); border = "replicate")

insp(0.0,1,0)


function f(center_val, periphery_val)
    center_val < periphery_val ? center_val : periphery_val
end


imgin = re[1:4, 1:4]
imgin[2,2] = 1.0f0
imgin .== extreme_filter(min, imgin, strel_box((3,3)))


insp(0,0,1)
result = imgin .== extreme_filter(f, imgin, ImageMorphology.StructuringElements.strel_box((3,3)))





# These are the sizes of (in pixels) of the component at this location
areas(MaxTree(re))
diameters(MaxTree(re))
# The maxima are numbered uniquely (except when neighbouring pixels have an identical value),
# but without a system we can use directly.
for i in eachindex(gr)
    if gr[i] > 0
        local_maxima_no = gr[i]
        local_maxima_value = re[i]
        print(Int(local_maxima_no), ':', Int(round(local_maxima_value, digits=0)), "   ")
    end
end
# The pixel areas associated with local maximaa are all 1.
matareas = areas(MaxTree(re))

# We can't see a direct way to extract saddle point (water shed) levels from the MaxTree. Saddle point levels
# should be those levels where the number of connected components in a max tree increases.
# Let's try to find candidates for inspection instead:
#
# There are lots of very close local maxima.
maxvals = sort([re[i] for i in eachindex(gr) if gr[i] > 0])
# We'll use the sorted maxvals to determine histogram bins.
saddles = maxvals .+ 0.5 * minimum(diff(maxvals))
edges = sort(vcat(maxvals, saddles))
# Let's have a look at the distribution of levels:
plo = zeros(Gray{Bool}, Int(round(maximum(edges) - minimum(edges) + 1)), length(edges))
for i in eachindex(edges)
    plo[Int(round(edges[i] - minimum(edges) + 1)), i] = Gray{Bool}(true)
end
plo
# Ok, there are just too many levels to consider (though considerably less than infinity).
# The largest level difference is
maximum(diff(maxvals))

# The procedure we have roughly in mind is:
# Use maxtree, reversed order.
# Assign a bool img for checking connectivity.
# Loop (maximano, z, parentindices) over corresponding arrays
#    Continue if maximano is zero
#    We can assume that the parent's elevation is larger.
#        els = filter 'edges' >= z
#        Loop through the lower elevations zl in els
#            Mask everyting below current mask level. This is a Boolean image M. Use in-place.
#            labeled_M, _ = label_components(M)
#            Is the local maximum connected to its parent? labeled_M[idx1] == labeled_M[idx16]
#                 Δz = local maxima - current mask level
#                 Mark Δz as the local maximum's prominence.

