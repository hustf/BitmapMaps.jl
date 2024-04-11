# Inspired by Luxor.Partition, adapted so as to ensure that an integer number of sheets
# cover the entire area, and that the oposite direction of 'northing' and 'image matrix y+'
# does not cause inaccuracies (which it might do). Another cause of bugs is the integer
# sampling distance.
#

"""
    SheetPartition

# Fields
- pixel_origin_ref_to_bitmapmap::Tuple{Int64, Int64}

The sheet's upper left corner given in a pixel coordinate system for the entire 'bitmapmap' (BM, consisting of multiple sheets).

- pix_iter::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}

Iterate over every pixel in the current sheet, starting top left. Typical use with `zip(pix_iter, utm_iter)`

- utm_iter::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}

Iterate over every utm coordinate (easting, northing) corresponding to a pixel in the current sheet, starting top left. Typical use with `zip(pix_iter, utm_iter)`


# Example

```
```
"""
@kwdef struct SheetPartition
    pixel_origin_ref_to_bitmapmap::Tuple{Int64, Int64}
    pix_iter::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}
    f_I_to_utm::Function
end

function SheetPartition(sheet_lower_left_utm, pix_dist, pixel_origin_ref_to_bitmapmap, pix_iter)
    sheet_height = size(pix_iter)[1]
    function f_I_to_utm(I::CartesianIndex)
        easting_offset =   (I[2] - 1) * pix_dist 
        northing_offset =  (sheet_height - I[1]) * pix_dist
        (easting_offset, northing_offset) .+ sheet_lower_left_utm
    end
    SheetPartition(;pixel_origin_ref_to_bitmapmap, pix_iter, f_I_to_utm)
end

function Base.show(io::IO, ::MIME"text/plain", shi::SheetPartition)
    println(io, "SheetPartition(")
    println(io, rpad("\tpixel_origin_ref_to_bitmapmap = ", 32), shi.pixel_origin_ref_to_bitmapmap)
    println(io, rpad("\tpix_iter or I = ", 32), shi.pix_iter)
    println(io, rpad("\tf_I_to_utm(first(pix_iter)) = ", 32), shi.f_I_to_utm(first(shi.pix_iter)), " )")
end


"""
Unless the entire 'bitmapmap' (the entire map or canvas) and 'sheet' 
dimensions are exact multiples of each other, the bordering sheets will  
extend outside of the 'bitmapmap' with background content.

TODO: Hence, we need to define a background color.
TODO: Input:
    1)  utm bounding box, pixels sheet width and height, and an integer utm_to_pixel_factor.
        Output is number of sheets.
    2) utm south-west, pixels sheet width and height, number of sheets m and n, and an integer utm_to_pixel_factor.

    Point is also used for utm coordinates elsewhere, though that is float.

    numbering: South-west to south-east, ending at north-east.
               This is nice, almost necessary because the source elevation grid is measured this way.
               Hence, we don't lose data at the north border.

    Output from each iteration:
    (m, n) bb_utm bb_pix

    pixel
"""
mutable struct BmPartition
    bm_pixel_width::Int
    bm_pixel_height::Int
    sheet_width::Int
    sheet_height::Int
    nrows::Int
    ncols::Int
    bm_southwest_corner::Tuple{Int, Int}
    sheet_indices::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}
    pixel_distance::Int
    function BmPartition(bm_pixel_width, bm_pixel_height, sheet_width, sheet_height, bm_southwest_corner, pixel_distance)
        ncols = bm_pixel_width / sheet_width |> ceil |> Integer
        nrows = bm_pixel_height / sheet_height |> ceil |> Integer
        if ncols < 1 || nrows < 1
            throw(error("partition(): not enough space for sheets that size"))
        end
        sheet_indices = CartesianIndices((1:nrows, 1:ncols))
        new(bm_pixel_width, bm_pixel_height, sheet_width, sheet_height, nrows, ncols, bm_southwest_corner, sheet_indices, pixel_distance)
    end
end

function SheetPartition(p::BmPartition, sheetnumber::Int)
    # No upper bounds here, 'iterate' and 'getindex' should do that.
    #
    # The first (and any) sheet's own pixel origin is top left of the sheet.
    #
    # The bitmapmap's pixel origin is top left of the bitmapmap (bm). 
    # which contains the south-west corner of the bitmapmap
    pix_iter = CartesianIndices((1:p.sheet_height, 1:p.sheet_width))
    r, c = row_col_of_sheet(p, sheetnumber)
    opx = (c - 1) * p.sheet_width
    opy = p.bm_pixel_height - p.sheet_height - (r - 1) * p.sheet_height
    pixel_origin_ref_to_bitmapmap = (opx, opy)
    sheet_lower_left_utm = p.bm_southwest_corner .+ p.pixel_distance .* (opx, (r - 1) * p.sheet_height)
    SheetPartition(sheet_lower_left_utm, p.pixel_distance, pixel_origin_ref_to_bitmapmap, pix_iter)
end

# The first iteration, returns first and second state.
function Base.iterate(p::BmPartition)
    sheetnumber = 1
    shp = SheetPartition(p, 1)
    shp_n = SheetPartition(p, 2)    
    return ((shp, sheetnumber), (shp_n, sheetnumber + 1))
end

# Other than the first iteration, returns this and the next state.
function Base.iterate(p::BmPartition, state)
    if state[2] > (p.nrows * p.ncols)
        return nothing
    end
    sheetnumber = state[2]
    shp = state[1]
    @assert shp isa SheetPartition
    shp_n = SheetPartition(p, sheetnumber + 1)
    return ((shp, sheetnumber), (shp_n, sheetnumber + 1))
end


function Base.length(p::BmPartition)
    p.nrows * p.ncols
end

row_col_of_sheet(p::BmPartition, sheetnumber::Int) = row_col_of_sheet(p.nrows, sheetnumber)
function row_col_of_sheet(nrows::Int, sheetnumber::Int)
    row = mod1(sheetnumber, nrows)
    col = div(sheetnumber - 1, nrows) + 1
    row, col
end

function Base.getindex(p::BmPartition, i::Int)
    1 <= i <= p.ncols * p.nrows || throw(BoundsError(p, i))
    shp = SheetPartition(p, i)
    return (shp, i)
end

#= For TODO
# TODO: Drop the 'ShIterators' name. Create a function returning 'i' from a ShIterator.
# See if we can't drop the 'i' from iterator returns.
function index_number_of_sheet(XXX::SheetPartition)
    2 
end
=#