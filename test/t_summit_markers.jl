using Test
using BitmapMaps
using BitmapMaps: CONSOLIDATED_FNAM, MaxTree, maximum_elevation_above, leaf_indices
using BitmapMaps: line!, mark_at!, prominence, distinct_summit_indices
using BitmapMaps: parent_of_leaf_indices, core_family_dictionary, get_consistent_random_color
using BitmapMaps: find_prominence
import ImageCore
using ImageCore: RGB, N0f8, RGBA, Gray
import Random

#################
# Unit test cases
#################
#  Case A
#
z_ = [1 1 1
      1 2 1
      1 1 1]
maxtree_ = MaxTree(z_)
parent_set_ = Set(maxtree_.parentindices)
@test parent_set_ == Set(1)
@test leaf_indices(maxtree_; parent_set = parent_set_) == 2:9
@test parent_of_leaf_indices(maxtree_; parent_set = parent_set_) == [1]
@test core_family_dictionary(maxtree_) == Dict(1 => [2, 3, 4, 5, 6, 7, 8, 9])
si_ = distinct_summit_indices(z_, maxtree_)
@test si_  == Set(5)
# We can't do this with Int64
@test_throws ArgumentError maximum_elevation_above(z_ ; maxtree = maxtree_, summit_indices = si_)
#
# Case B
# 
z_ =  Float32[ 1.0  1.0  1.0
        1.0  2.0  2.2
        1.0  2.1  1.0]
maxtree_ = MaxTree(round.(z_))
parent_set_ = Set(maxtree_.parentindices)
@test parent_set_ == Set([1, 5])
@test leaf_indices(maxtree_; parent_set = parent_set_) == [2, 3, 4, 6, 7, 8, 9]
@test parent_of_leaf_indices(maxtree_; parent_set = parent_set_) == [1, 5]
@test core_family_dictionary(maxtree_) == Dict(5 => [6, 8])
si_ = distinct_summit_indices(z_, maxtree_)
@test si_  == Set(8)
@test ! any(isnan, maximum_elevation_above(z_ ; maxtree = maxtree_, summit_indices = si_))

#
# Case 
# 
z_ = Float32[ 21.0  19.0  19.0  16.0
       22.0  22.0  23.0  22.0
       22.0  23.0  23.0  21.0
       22.0  22.0  22.0  16.0]
maxtree_ = MaxTree(round.(z_))
parent_set_ = Set(maxtree_.parentindices)
@test parent_set_ == Set([1, 2, 5, 7, 13])
@test core_family_dictionary(maxtree_) == Dict(7 => [10, 11])
si_ = distinct_summit_indices(z_, maxtree_)
@test si_   == Set(7)
@test ! any(isnan, maximum_elevation_above(z_ ; maxtree = maxtree_, summit_indices = si_))

#
# Case 
# 
z_ =  Float32[ 1.0  1.0  1.0  2.2
        1.0  2.0  1.2  3.0
        1.0  2.1  1.0  1.0]
maxtree_ = MaxTree(round.(z_))
parent_set_ = Set(maxtree_.parentindices)
@test parent_set_ == Set([1, 5, 10])
@test sort(leaf_indices(maxtree_; parent_set = parent_set_)) == [2, 3, 4, 6, 7, 8, 9, 11, 12]
@test parent_of_leaf_indices(maxtree_; parent_set = parent_set_) == [1, 5, 10]
@test core_family_dictionary(maxtree_) == Dict(5 => [6], 10 => [11]) 
si_ = distinct_summit_indices(z_, maxtree_)
@test si_  == Set([6, 11])
@test ! any(isnan, maximum_elevation_above(z_ ; maxtree = maxtree_, summit_indices = si_))

# Case 
# Peaks: 7 and 15
z_ = Float32[  0.0   1.0   0.0  -3.0   -8.0
                1.0   2.0   1.0   2.5   -7.0
                0.0   2.1   0.0  10.0   -8.0
              -3.0  -2.0  -3.0  -6.0  -11.0]
maxtree_ = MaxTree(round.(z_))
@test length(filter(i -> maxtree_.parentindices[i] == i, maxtree_.traverse)) == 1
parent_set_ = Set(maxtree_.parentindices)
@test core_family_dictionary(maxtree_) == Dict(6 => [7], 14 => [15])
si_ = distinct_summit_indices(z_, maxtree_)
@test si_  == Set([7, 15])
@test ! any(isnan, maximum_elevation_above(z_ ; maxtree = maxtree_, summit_indices = si_))


# Case 
# Peaks: 6 and 15
z_ = Float32[  0.0   1.0   0.0  -3.0   -8.0
               1.0   2.0   1.0  -2.0   -7.0
               0.0   1.0   0.0   4.0    3.0
              -3.0  -2.0  -3.0  -6.0  -11.0]
maxtree_ = MaxTree(round.(z_))
parent_set_ = Set(maxtree_.parentindices)
@test core_family_dictionary(maxtree_) == Dict(2 => [5, 6, 7, 10], 19 => [15])
si_ = distinct_summit_indices(z_, maxtree_)
@test si_  == Set([6, 15])
@test ! any(isnan, maximum_elevation_above(z_ ; maxtree = maxtree_, summit_indices = si_))



################
# Assembly tests
################
fofo = joinpath(homedir(), "BitmapMaps", "proj 47675 6929520 57224 6947852", "1 1  47675 6929520  57224 6943269")
if isdir(fofo)  # Created with cell2utm = 3
#
@test isfile(joinpath(fofo, CONSOLIDATED_FNAM))

# Read the consolidated file and confirm the utm location function
g = readclose(joinpath(fofo, CONSOLIDATED_FNAM))
utmne = Int.(northeast_external_corner(g))
@test parse.(Int, split(splitdir(fofo)[end])[5:6]) .- utmne == [0, 0]
utmsw = Int.(southwest_external_corner(g))
@test parse.(Int, split(splitdir(fofo)[end])[3:4]) .- utmsw == [0, 0]
utmcen = (utmne .+ utmsw) .÷ 2
icen_g = length(g.A) ÷ 2
@test Int.(g.f(CartesianIndices(g.A[:, :, 1])[icen_g].I)) .- utmcen == [0, 0]

# Downsample elevations z
cell_iter = CartesianIndices((1:4583, 1:3183))
cell2utm = 3
ny, nx = size(cell_iter)
source_indices = (1:cell2utm:(nx  * cell2utm), 1:cell2utm:(ny * cell2utm))
si = CartesianIndices(source_indices)
z = transpose(g.A[si])

# Find and check the maxtree
maxtree = MaxTree(round.(z))
@test length(filter(i -> maxtree.parentindices[i] == i, maxtree.traverse)) == 1
# Find the tallest index from each summit.
summit_indices = distinct_summit_indices(z, maxtree)
# Find and check maximum elevation above
mea = maximum_elevation_above(z; maxtree, summit_indices)
@test ! any(isnan, mea)
# Find and check prominence. Since we don't consider neighbouring sheets,
# and the minimum elevation in this sheet is zero, the prominence of the tallest peak equates
# its elevation.
prom = prominence(z, summit_indices, mea; maxtree)
@test maximum(filter(p -> ! isnan(p), prom)) == 1431.9824f0
# Inspect height and mark maximum
display_if_vscode(z)
maxel, maxind = findmax(z)
mark_at!(z, [maxind], 601, "on_circle")
display_if_vscode(z)
# Get rid of the marker
z = transpose(g.A[si])

#
# The following ought to run both with and without pre-existing file '_Max_elevation_above.mat'
# Manual deletion is OK. 
#

# Test and inspect "Randers topp"
I = CartesianIndex(1352, 1500)
@test abs(z[I] - 1414.5989f0) < 0.01 
@test abs(prom[I] - 137.00366f0) < 0.01
@test isnan(prom[I + CartesianIndex(1,1)])
hw = 5
I_close = CartesianIndices((I[1] - hw:I[1] + hw, I[2] - hw : I[2] + hw))
display_if_vscode(z[I_close])

# Test and inspect "Ramoen"
I = CartesianIndex(1227 - 5 + 0, 1369 + 5 + 2 + 1)
I_lin = LinearIndices(z)[I]
@test I_lin ∈ summit_indices
@test abs(z[I] - 1418.0859f0) < 0.01 
@test abs(prom[I] - 453.72192f0) < 0.01
@test isnan(prom[I + CartesianIndex(1,1)])
hw = 5
I_close = CartesianIndices((I[1] - hw:I[1] + hw, I[2] - hw : I[2] + hw))
display_if_vscode(z[I_close])


##########################################
# Distinguish clearly between summit zones
# by inspecting indexed image
##########################################
imea = Int.(round.(mea))
levcols = get_consistent_random_color.(minimum(imea) : maximum(imea))

indimg = levcols[imea]
hw = 500
I_close = CartesianIndices((I[1] - hw:I[1] + hw, I[2] - hw : I[2] + hw))
display_if_vscode(indimg[I_close])
# Save for inspection
ffna = joinpath(fofo, "summit_regions.png")
save_png_with_phys(ffna, indimg)

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


##############################################
# Trace such a path for all prominent summits.
# Traces illustrate the maxtree.
##############################################

traces = zeros(Float32, size(z)...)
inds = map(i -> CartesianIndices(z)[i], findall(p -> ! isnan(p) && p > 100, prom))
filter!(i -> z[i] > 200, inds)
@test length(inds) < 50
# For efficiency, sort inds by elevation - tallest first.
vz = z[inds]
order = sortperm(vz; rev = true)
mark_at!(traces, inds[order]; side = 49, f_is_filled = BitmapMaps.func_is_in_triangle)
display_if_vscode(traces)
fill!(mea, NaN)
for ind in inds[order]
    el = z[ind]
    parent_i = maxtree.parentindices[ind]
    print(".") # Indicate progress. The first dots are the slowest to appear, faster towards end.
    trace_downpath!(traces, mea, maxtree, el, parent_i, ind)
end
display_if_vscode(traces)

# Save for inspection
ffna = joinpath(fofo, "traces_downhill.png")
save_png_with_phys(ffna, map(traces) do pix
    pix == true && return RGBA{N0f8}(0, 0, 0, 1)
    RGBA{N0f8}(0., 0, 0, 0)
end)





# Filter out obscure summits that are really artifacts or
# power line spans.

g11, _, _, g22 = (BitmapMaps.hessian_components(collect(z)))
indicator = map(g11 .+ g22) do σ
    Float32(σ < 0 ? σ : 0.0f0)
end

indicator1 = map(g11 .+ g22) do σ
    Float32(σ)
end
display_if_vscode(indicator ./ maximum(indicator))
minimum(indicator)
# Save for inspection
ffna = joinpath(fofo, "principal_tensile_stress.png")
save_png_with_phys(ffna, map(indicator ./ minimum(indicator)) do pix
    Gray{N0f8}(pix)
end)


prom[cell_iter]
z[cell_iter]

peak_stress = [prom[i] > 30 ? indicator1[i] : 0.0f0 for i in cell_iter]
summit_elevation = [prom[i] > 30 ? z[i] : 0.0f0 for i in cell_iter]
minimum(peak_stress)
maximum(peak_stress)
sort(filter(s -> s < 0 ,peak_stress))


else
    @info "Skipping assembly tests because data is missing"
end # isdir, after_assembly_test
