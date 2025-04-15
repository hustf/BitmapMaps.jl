# Utiltiy functions for GeoArrays.
# Some extend functions defined for SheetBuilder and SheetMatrixBuilder


"""
    readclose(fna)

A barrier for a possible issue with GeoArrays. Keyword arguments not included.
Because calling `read` sometimes leaves a file 'open', runs in a separate thread.
"""
function readclose(fna)::GeoArray # Type annotation for the compiler, just like `read`
    file_task = Threads.@spawn GeoArrays.read(fna)
    contents = fetch(file_task)
    GC.gc()            # Clean up resources explicitly after the operation
    # We can't make stirct type assertions here, because some unconsolidated files return
    # GeoArrays.GeoArray{Union{Missing, Float32}, Array{Union{Missing, Float32}, 3}}
    # Do make type assertion in the caller.
    # @assert contents isa GeoArray{Float32, Array{Float32, 3}} "$(typeof(contents))"
    contents
end

# See docstring in builders_utilties, this extends.
function closed_polygon_string(fnas::Vector{String})
    s = ""
    n = length(fnas)
    # Let's avoid taking several g into memory simultaneously
    for (i, fna) in enumerate(fnas)
        @debug "Finding extents of $fna"
        g = readclose(fna)
        @assert g isa GeoArray{Float32, Array{Float32, 3}} "$(typeof(g))"
        s *= closed_polygon_string(g)
        if i < n
            s *= ",\n" * repeat(' ', 19)
        end
    end
    s
end
function closed_polygon_string(fna::String)
    # This does not normally occur in the pipeline, rather while preparing input.
    @debug "Finding extents of $fna"
    closed_polygon_string(readclose(fna))
end
function closed_polygon_string(g::GeoArray)
    bbo = bbox(g)
    bbi = nonzero_raster_rect(g)
    so = "$(closed_box_string(bbo))"
    if bbi == (min_x = 0, min_y = 0, max_x = 0, max_y = 0)
        # This file is all zero-valued or missing-valued.
        # We indicate that visually with a diagonal line across the outer
        # rectangle.
        si = "$(diagonal_string(bbo))"
        so * ",\n" * repeat(' ', 19) * si
    elseif bbi == bbo
        # This file has non-zero data filling the outer boundaries.
        so
    else
        # The inner rectangle will indicate where non-zero data occurs.
        si = "$(closed_box_string(bbi))"
        so * ",\n" * repeat(' ', 19) * si
    end
end

"""
    _convToNamedTuple(bb::Extents.Extent)

We rely on GeoArrays.jl, which recently was revised to use Extents.jl instead of
NamedTuple. This conversion is part of a minimal change solution.
"""
function _convToNamedTuple(bb::Extents.Extent)
    #T = NamedTuple{(:min_x, :min_y, :max_x, :max_y)}
    #(min_x = bb.X[1], min_y = bb.Y[1], max_x = bb.X[2], max_y = bb.Y[2])
    NamedTuple{(:min_x, :min_y, :max_x, :max_y)}((bb.X[1], bb.Y[1], bb.X[2], bb.Y[2]))
end
# Docstring in builders_utilties
function bbox_external_string(fna::String)
    g = readclose(fna)
    @assert g isa GeoArray{Float32, Array{Float32, 3}} "$(typeof(g))"
    min_x, min_y = southwest_external_corner(g)
    max_x, max_y = northeast_external_corner(g)
    bbox_external_string((;min_x, min_y, max_x, max_y))
end

# Docstring in builders_utilties
cell_to_utm_factor(fna::String) = cell_to_utm_factor(readclose(fna))
function cell_to_utm_factor(g::GeoArray)
    g.f((1, 0))[1] - g.f((0, 0))[1]
end



"""
    nonzero_raster_rect(fna::String)
    nonzero_raster_rect(g::GeoArray)
    ---> @NamedTuple{min_x::Int64, min_y::Int64, max_x::Int64, max_y::Int64}

This returns the geographical external boundaries of the cells containing data other than zero.
It also stores the result in memory. See TIFDIC constant note.

I.e. unpadded external utm boundaries. See figure `/resource/padded_geoarray.svg`.
"""
function nonzero_raster_rect(fna::String)
    @assert isfile(fna)
    @assert endswith(fna, r".tif|.TIF")
    if haskey(TIFDIC, fna)
        bb = TIFDIC[fna]
    else
        bb = nonzero_raster_rect(readclose(fna))
        push!(TIFDIC, fna => bb)
    end
    bb
end
function nonzero_raster_rect(g::GeoArray)::@NamedTuple{min_x::Int64, min_y::Int64, max_x::Int64, max_y::Int64}
    # Downloaded .tifs are often zero padded.
    # Indices of outer cells that contain non-zero data.
    itup = unpadded_indices(g.A[:,:,1])
    if itup == (nothing, nothing)
        return (; min_x = 0, min_y =  0, max_x = 0, max_y = 0)
    end
    Irng = CartesianIndices(itup)
    x4, y4 = northwest_corner(g, Irng)
    x2, y2 = southeast_external_corner(g, Irng)
    x3, y3 = x2, y4 # NE
    x1, y1 = x4, y2 # SW
    # External' corner, or North-West corner of the NW cell.
    max_x, max_y = x2, y4
    min_x, min_y = x4, y2
    # return NamedTuple
    (; min_x, min_y, max_x, max_y)
end



"""
    nonzero_raster_closed_polygon_string(fna::String)

Unpadded (non-zero) exteral extents. See `closed_polygon_string` for builders.
"""
function nonzero_raster_closed_polygon_string(fna::String)
    bbup = nonzero_raster_rect(fna)
    closed_box_string(bbup)
end

"""
unpadded_indices(matrix)

When you want to drop all zero padding.

# Example
```
julia> M = [
    0 0 0 0
    0 1 0 0
    0 3 0 0
    0 0 0 0];

julia> Irng = BitmapMaps.unpadded_indices(M)
    CartesianIndices((2:3, 2:2))

julia> M[Irng]
    2×1 Matrix{Int64}:
     1
     3
```
"""
function unpadded_indices(matrix)
    # Get the number of rows and columns
    m, n = size(matrix)
    # Determine the first and last non-zero row
    first_row = findfirst(!is_all_zeros, eachrow(matrix))
    last_row = findlast(!is_all_zeros, eachrow(matrix))
    # Determine the first and last non-zero column
    first_col = findfirst(!is_all_zeros, eachcol(matrix))
    last_col = findlast(!is_all_zeros, eachcol(matrix))
    # Return the ranges for rows and columns to keep
    if isnothing(first_row) && isnothing(first_col)
        return(nothing, nothing)
    else
        (first_row:last_row, first_col:last_col)
    end
end
function unpadded_indices(g::GeoArray)
    is, js = unpadded_indices(g.A[:, :])
    @assert is !== nothing
    @assert js !== nothing
    CartesianIndices((is, js))
end

# Function to check if a row or column is all zeros
is_all_zeros(vec::SubArray{Float32}) = all(x -> x == 0, vec)
is_all_zeros(vec::SubArray{T}) where T <: Union{Missing, Float32} = all(x -> ismissing(x) || x == 0, vec)

"ref. southwest_corner"
northwest_corner(g::GeoArray) = g.f((0, 0))

function northwest_corner(g::GeoArray, Irng)
    g.f(Tuple(Irng[1]) .-(1,1))
end

"ref. southwest_corner"
southeast_internal_corner(g::GeoArray) = g.f(size(g)[1:2] .-(1,1))
southeast_internal_corner(g::GeoArray, Irng) = g.f(Tuple(Irng[end]) .-(1,1))

"ref. southwest_corner"
southeast_external_corner(g::GeoArray) = g.f(size(g)[1:2])
southeast_external_corner(g::GeoArray, Irng) = g.f(Tuple(Irng[end]))