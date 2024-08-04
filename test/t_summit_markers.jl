using Test
using BitmapMaps
using BitmapMaps: CONSOLIDATED_FNAM, MaxTree, maximum_elevation_above, leaf_indices
using BitmapMaps: line!, mark_at!, prominence_clusters
import ImageCore
using ImageCore: RGB, N0f8, RGBA
using ImageSegmentation
import Random
import ImageMorphology
using ImageMorphology: local_maxima

fofo = joinpath(homedir(), "BitmapMaps", "proj 47675 6929520 57224 6947852", "1 1  47675 6929520  57224 6943269")
@assert ispath(fofo)
@assert isfile(joinpath(fofo, CONSOLIDATED_FNAM))

# Read the consolidated file and confirm the utm location function
g = readclose(joinpath(fofo, CONSOLIDATED_FNAM))
utmne = Int.(northeast_external_corner(g))
@test parse.(Int, split(splitdir(fofo)[end])[5:6]) .- utmne == [0, 0]
utmsw = Int.(southwest_external_corner(g))
@test parse.(Int, split(splitdir(fofo)[end])[3:4]) .- utmsw == [0, 0]
utmcen = (utmne .+ utmsw) .÷ 2
icen_g = length(g.A) ÷ 2
@test Int.(g.f(CartesianIndices(g.A[:, :, 1])[icen_g].I)) .- utmcen == [0, 0]

# Downsample and round off elevations to whole metres.
cell_iter = CartesianIndices((1:4583, 1:3183))
cell2utm = 3
ny, nx = size(cell_iter)
source_indices = (1:cell2utm:(nx  * cell2utm), 1:cell2utm:(ny * cell2utm))
si = CartesianIndices(source_indices)
z = transpose(round.(g.A[si]))

# Extract maximum elevation above and prominence
maxtree = MaxTree(z)
mea = maximum_elevation_above(z; maxtree)
prom = prominence_clusters(z;maxtree, mea)

# Inspect height and mark maximum
display_if_vscode(z)
maxel, maxind = findmax(z)
mark_at!(z, [maxind], 601, "on_circle")
display_if_vscode(z)
# Get rid of the marker
z = transpose(round.(g.A[si]))


##########################################
# Distinguish clearly between summit zones
# by inspecting segmented image
##########################################
levels = sort(unique(mea))
seeds = map(lev -> (findfirst(z-> z==lev , mea), Int(lev)), levels)
seg = seeded_region_growing(mea , seeds)
function get_random_color(seed)
    Random.seed!(seed) # For consistentency between runs
    rand(RGB{N0f8})
end
segcol = map(i -> get_random_color(i), labels_map(seg));
display_if_vscode(segcol)
display_if_vscode(segcol[200:700, 1:500])
display_if_vscode(mea)
display_if_vscode(mea[200:700, 1:500])

###################################################
# Trace a path followed by algo max_elevation_above
###################################################

# This is a modification of function _mea!
function trace_downpath!(traces, mea, maxtree, leaf_elevation, i, child_i)
    previous_maximum_elevation_above = mea[i]
    if isnan(previous_maximum_elevation_above)
        # This line is just for understanding the algorithm.
        line!(traces, CartesianIndices(mea)[child_i], CartesianIndices(mea)[i])
        # Remember this leaf elevation. Most likely, another and taller
        # leaf will overwrite this later.
        mea[i] = leaf_elevation
    elseif previous_maximum_elevation_above >= leaf_elevation
        # This line is just for understanding the algorithm.
        line!(traces, CartesianIndices(mea)[child_i], CartesianIndices(mea)[i])
        # Just leave quietly. Stop the recursion here.
        return mea
    else
        # This line is just for understanding the algorithm.
        line!(traces, CartesianIndices(mea)[child_i], CartesianIndices(mea)[i])
        # Forget about what we thought was the summit elevation.
        # Remember this instead!
        mea[i] = leaf_elevation
    end
    # Let's look for our parent.
    parent_i = maxtree.parentindices[i]
    if parent_i == i
        # ⬆ We are our own parent, the root, lowest level in the map.
        # Stop the recursion here.
        return mea
    end
    # ⬇ Recurse: Tell current parent (which is lower down) 
    # about the leaf's elevation. It's the actual summit above us for all we know.
    trace_downpath!(traces, mea, maxtree, leaf_elevation, parent_i, i)
end


###########################################
# Trace such path for all prominent summits
###########################################

traces = zeros(Float32, size(z)...)
inds = map(i -> CartesianIndices(z)[i], findall(p -> ! isnan(p) && p > 100, prom))
mark_at!(traces, inds; side = 19, f_is_filled = BitmapMaps.func_is_in_triangle)
display_if_vscode(traces)
fill!(mea, NaN)
for ind in inds
    el = z[ind]
    parent_i = maxtree.parentindices[ind]
    trace_downpath!(traces, mea, maxtree, el, parent_i, ind)
    print(".")
end
display_if_vscode(traces)

# Save for inspection
ffna = joinpath(fofo, "traces_downhill.png")
save_png_with_phys(ffna, map(traces) do pix
    pix == true && return RGBA{N0f8}(0, 0, 0, 1)
    RGBA{N0f8}(0., 0, 0, 0)
end)
