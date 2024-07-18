# NOTE:
# This is work in progress. Also related: t_peaks.jl and mark_utils.jl. We are trying to identify and mark 
# peaks from their prominence ("primary factor"), and while to reusing some functions from 'contour.jl'.
# Marking may be delegated to the svg.
using Test
using BitmapMaps
using BitmapMaps: indices, coords
import ColorTypes
using ColorTypes: RGB, N0f8, Gray
using ImageView
using ImageCore: channelview, zeroarray, colorview, scaleminmax
import ImageMorphology
using ImageMorphology: MaxTree, local_maxima, label_components, areas, label_components!
#import ImageSegmentation
#using ImageSegmentation: 

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
insp(1,1,1)

# The procedure we have roughly in mind is:
# Use maxtree, reversed order.
# Assign a bool img for checking connectivity.
# Loop (max_index, z, parent) over corresponding arrays
#    Continue if maximano is zero
#    We can assume that the parent's elevation is larger.
#        elevs = filter 'edges' >= z
#        Loop through the lower elevations zl in elevs 
#            Mask everyting below current mask level. This is a Boolean image M. Use in-place.
#            labeled_M, _ = label_components(M)
#            Is the local maximum connected to its parent? labeled_M[idx1] == labeled_M[idx16]
#                 Δz = local maxima - current mask level
#                 Mark Δz as the local maximum's prominence.

mtree = MaxTree(re, rev = true)
loc_ma = local_maxima(re)
function candidate_elevations(loc_ma, z_ma)
    # NOTE this is flawed. Maybe we ought to use local minima here instead.
    # Or we could identify all saddle point elevations, see ridges.jl.
    # Those could be used directly.
    @assert size(loc_ma) == size(z_ma)
    maxvals = sort([z_ma[i] for i in eachindex(loc_ma) if loc_ma[i] > 0])
    Δ = 0.5f0 * minimum(diff(maxvals))
    Δ <= 0 && throw("Unexpected non-unique maxvals. Call unique first if this happens.")
    sort(vcat(maxvals, maxvals .+ Δ), rev = true)
end
elevs = candidate_elevations(loc_ma, re)
#=
bl .= scaleminmax(extrema(areas(mtree))...).(areas(mtree))
insp(1,0,1)
label_components(M)
gr .= scaleminmax(extrema(areas(mtree))...).(areas(mtree))
=#

bwbuf = zeros(Gray{Bool}, size(bl))
prominence = zeros(Gray{Float32}, size(bl))
labelbuf = zeros(Int64, size(bl))
for (locmax_i, z, parent_i) in zip(eachindex(loc_ma), re, mtree.parentindices)
    iszero(locmax_i) && continue
    for zlevel in elevs
        zlevel >= z && continue
        # In bwbuf, everything below zlevel is false.
        # Everything else is an island.
        # This creates more or less isolated islands. 
        map!(zl -> Gray{Bool}(zl > zlevel), bwbuf, re)
        bwbuf .= Gray{Bool}.(re .> zlevel)
        # Check if zlevel will connect locmax_i to parent_i
        label_components!(labelbuf, bwbuf) # TODO in-place!
        #@show locmax_i parent_i
        #display(bwbuf)
        if labelbuf[locmax_i] == labelbuf[parent_i]
            # This is the definition of prominence: How much lower did we 
            # go to connect to a higher level?
            Δz = z - zlevel
            if Δz >= 30
                prominence[locmax_i] = Δz
                # TODO: Make a triangular 'false' mark in bwbuf 
                # Make a square 'false' mark in bw_buf
                #coords(g, (nw_ind[1]  + ii, nw_ind[2] + jj))
                @show z Δz 
                display(transpose(bwbuf))
            end
            # We're through - finish the elevation loop and contnue with the next local maximum
            break
        end
        if zlevel == last(elevs)
            @show "Found no prominence for local maximum at z = $z"
        end
    end
end