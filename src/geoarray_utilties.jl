"""
    readclose(fnam)

A barrier for a possible issue with GeoArrays. Keyword arguments not included.
 Because calling `read` sometimes leaves a file 'open', runs in a separate thread.
"""
function readclose(fnam)::GeoArray # Type annotation for the compiler, just like `read` 
    file_task = Threads.@spawn GeoArrays.read(fnam)
    contents = fetch(file_task)
    GC.gc()            # Clean up resources explicitly after the operation
    contents
end



"""
    nonzero_raster_rect(fna::String)::@NamedTuple{min_x::Int64, min_y::Int64, max_x::Int64, max_y::Int64}

This returns the geographical external boundaries of the cells containing data other than zero.

I.e. unpadded external utm boundaries. See figure `/resource/padded_geoarray.svg`.
"""
function nonzero_raster_rect(fna::String)
    @assert isfile(fna)
    @assert endswith(fna, r".tif|.TIF")
    nonzero_raster_rect(readclose(fna))
end
function nonzero_raster_rect(g::GeoArray)::@NamedTuple{min_x::Int64, min_y::Int64, max_x::Int64, max_y::Int64}
    # Downloaded .tifs are often zero padded.
    # Indices of outer cells that contain non-zero data.
    Irng = CartesianIndices(unpadded_indices(g.A[:,:,1]))
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
    # Function to check if a row or column is all zeros
    is_all_zeros = vec -> all(x -> x == 0, vec)
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
    is, js = unpadded_indices(g.A[:, :, 1])
    @assert is !== nothing
    @assert js !== nothing
    CartesianIndices((is, js))
end

"ref. southwest_corner"
northwest_corner(g::GeoArray) = g.f((1, 1))
function northwest_corner(g::GeoArray, Irng) 
    g.f(Tuple(Irng[1]))
end
"ref. southwest_corner"
southeast_internal_corner(g::GeoArray) = g.f(size(g)[1:2])
southeast_internal_corner(g::GeoArray, Irng) = g.f(Tuple(Irng[end]))
"ref. southwest_corner"
southeast_external_corner(g::GeoArray) = 2 .* southeast_internal_corner(g) .- g.f(size(g)[1:2] .-(1,1))
southeast_external_corner(g::GeoArray, Irng) = 2 .* southeast_internal_corner(g, Irng) .- g.f(Tuple(Irng[end]) .- (1, 1))
