using Test
using BitmapMaps
using BitmapMaps: CONSOLIDATED_FNAM, MaxTree, maximum_elevation_above, leaf_indices
using BitmapMaps: line!, mark_at!, prominence
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

# Check the maxtree
maxtree = MaxTree(z)
@test length(filter(i -> maxtree.parentindices[i] == i, maxtree.traverse)) == 1
# Extract maximum elevation above and prominence
mea = maximum_elevation_above(z; maxtree)
prom = prominence(z;maxtree, mea)

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

# Investigate why this occurs:
# Prominence_m        	Elevation_m         	Utm                 	Sheet_indices       
# 1432.0              	1432.0              	(49742, 6933327)    	(3315, 690)         
# 1432.0              	1432.0              	(49745, 6933312)        (3320, 691)   
I1 = CartesianIndex(3315, 690)
I2 = CartesianIndex(3320, 691)
@test z[I1] == z[I2]


# Check that both are leaf nodes:
parent_set = Set(maxtree.parentindices)
@test I1 ∉ parent_set
@test I2 ∉ parent_set
# Are they both children of the same node?
parent_i1 = maxtree.parentindices[I1]
parent_i2 = maxtree.parentindices[I2]
@test parent_i1 == parent_i2
# This is key! We need to refine the selection method.
# Just pick the tallest leaf among these....

li = leaf_indices(maxtree)
@test length(li) == length(unique(li))
leaf_parents = maxtree.parentindices[li]
counts = Dict{Int, Int}()
# Count occurrences of each parent
for x in leaf_parents
    counts[x] = get(counts, x, 0) + 1
end
# Parents that sire multiple leafs
repeat_parents = [k for (k, v) in counts if v > 1]
# # Parents with only one leaf
single_leaf_parents = setdiff(leaf_parents, repeat_parents)
@test length(single_leaf_parents) + length(repeat_parents) < length(leaf_parents)