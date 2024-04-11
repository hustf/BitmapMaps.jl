#=
Inspired by Luxor.Partition, adapted so as to ensure that 
 - an integer number of sheets covers the entire area
 - the oposite direction of 'northing' and 'image matrix y+' does not cause inaccuracies (which it might do). 
 - Hide away the bug nest of sampling distances.
=#

"""
    SheetPartition

A SheetPartition is one sheet in a collection defined by `BmPartition`. 

# Fields
- pixel_origin_ref_to_bitmapmap::Tuple{Int64, Int64}

The sheet's upper left corner given in a pixel coordinate system for the entire 'bitmapmap' (BM, consisting of multiple sheets).

- pix_iter::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}

Iterate over every pixel index `I` in the current sheet, starting top left. This field is the same for every sheet in a BmPartition.

- f_I_to_utm::Function  

f_I_to_utm(I) transforms a pixel index I for this sheet to an UTM coordinate. For UTM coordinates, see `GeoArrays.jl`.

- sheet_number::Int

Sequential index `i` of this sheet. You can get the sheet's row and column by calling `row_col_of_sheet(bmp, i)` 


# Example

```
julia> using BitmapMaps

julia> bitmapmap = let
           bm_pixel_width, bm_pixel_height = 9, 8  # Actual full bitmapmap sizes would be thousands of pixels
           sheet_width, sheet_height = 3, 4        # 12 sheets, over which bm_pixels can be evenly divided
           bm_southwest_corner = (43999, 6909048)  # Utm zone coordinates, the corner is lower left in bmp[1, 1] (defined below)
           pixel_distance = 2                      # One horizontal pixel distance equals 2 metres easting and 2 metres northing
           #
           BmPartition(bm_pixel_width, bm_pixel_height, sheet_width, sheet_height, bm_southwest_corner, pixel_distance)
       end;


julia> display.(bitmapmap );   # Broadcasting, mapping, other iterations are defined.
       SheetPartition(
               pixel_origin_ref_to_bitmapmap = (0, 4)
               pix_iter or I =                 CartesianIndices((1:4, 1:3))
               f_I_to_utm(first(pix_iter)) =   (43999, 6909054)
               sheet_number =                  1
        )
        ⋮
               sheet_number =                  12
        )

julia> bitmapmap[6] |> println     # Single index, short form printing.
SheetPartition((6, 0), (1:4, 1:3), f(I) -> utm, 6 )

julia> bitmapmap[2, 3] |> println  # Row, col index, same object
SheetPartition((6, 0), (1:4, 1:3), f(I) -> utm, 6 )

julia> utm_coordinates_of_pixels_in_sheet = let
    sheetpartition = bitmapmap[1]
    map(sheetpartition.f_I_to_utm, sheetpartition.pix_iter)
end
4×3 Matrix{Tuple{Int64, Int64}}:
(43999, 6909054)  (44001, 6909054)  (44003, 6909054)
(43999, 6909052)  (44001, 6909052)  (44003, 6909052)
(43999, 6909050)  (44001, 6909050)  (44003, 6909050)
(43999, 6909048)  (44001, 6909048)  (44003, 6909048)
```
"""
@kwdef struct SheetPartition
    pixel_origin_ref_to_bitmapmap::Tuple{Int64, Int64}
    pix_iter::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}
    f_I_to_utm::Function
    sheet_number::Int
end

function SheetPartition(sheet_lower_left_utm, pix_dist, pixel_origin_ref_to_bitmapmap, pix_iter, sheet_number)
    sheet_height = size(pix_iter)[1]
    function f_I_to_utm(I::CartesianIndex)
        easting_offset =   (I[2] - 1) * pix_dist 
        northing_offset =  (sheet_height - I[1]) * pix_dist
        (easting_offset, northing_offset) .+ sheet_lower_left_utm
    end
    SheetPartition(;pixel_origin_ref_to_bitmapmap, pix_iter, f_I_to_utm, sheet_number)
end

function Base.show(io::IO, ::MIME"text/plain", shi::SheetPartition)
    colwi = 32
    println(io, "SheetPartition(")
    println(io, rpad("\tpixel_origin_ref_to_bitmapmap = ", colwi), shi.pixel_origin_ref_to_bitmapmap)
    println(io, rpad("\tpix_iter or I = ", colwi), shi.pix_iter)
    println(io, rpad("\tf_I_to_utm(first(pix_iter)) = ", colwi), shi.f_I_to_utm(first(shi.pix_iter)))
    println(io, rpad("\tsheet_number = ", colwi), shi.sheet_number)
    println(io, " )")
end
function Base.show(io::IO, shi::SheetPartition)
    print(io, "SheetPartition(")
    print(io, shi.pixel_origin_ref_to_bitmapmap, ", ")
    print(io, shi.pix_iter.indices, ", ")
    print(io, "f(I) -> utm", ", ")
    print(io, shi.sheet_number)
    println(io, " )")
end



"""
    BmPartition

# Fields



See `SheetPartition` and `resource/map_sheet_utm_pix` for an example.
"""
struct BmPartition
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


function Base.show(io::IO, ::MIME"text/plain", p::BmPartition)
    colwi = 32
    println(io, "BmPartition(")
    for sy in [:bm_pixel_width,
            :bm_pixel_height,
            :sheet_width,
            :sheet_height,
            :nrows,
            :ncols,
            :bm_southwest_corner,
            :sheet_indices,
            :pixel_distance]
        println(io, rpad("\t$sy", colwi), " = ", getfield(p, sy), ",")
    end
    println(io, " )")
end


function _SheetPartition(p::BmPartition, sheet_number::Int)
    # No bounds in this "private" function, 'iterate' and 'getindex' should do that.
    pix_iter = CartesianIndices((1:p.sheet_height, 1:p.sheet_width))
    r, c = row_col_of_sheet(p, sheet_number)
    opx = (c - 1) * p.sheet_width
    opy = p.bm_pixel_height - p.sheet_height - (r - 1) * p.sheet_height
    pixel_origin_ref_to_bitmapmap = (opx, opy)
    sheet_lower_left_utm = p.bm_southwest_corner .+ p.pixel_distance .* (opx, (r - 1) * p.sheet_height)
    SheetPartition(sheet_lower_left_utm, p.pixel_distance, pixel_origin_ref_to_bitmapmap, pix_iter, sheet_number)
end

Base.iterate(bmp::BmPartition) = _SheetPartition(bmp, 1), _SheetPartition(bmp, 2)

function Base.iterate(p::BmPartition, state::SheetPartition)
    state.sheet_number > p.nrows * p.ncols && return nothing
    # return: the one to use next, the one after that
    state, _SheetPartition(p, state.sheet_number + 1)
end


Base.length(bmp::BmPartition) = bmp.nrows * bmp.ncols

row_col_of_sheet(p::BmPartition, sheet_number::Int) = row_col_of_sheet(p.nrows, sheet_number)
row_col_of_sheet(nrows::Int, sheetnumber::Int) = mod1(sheetnumber, nrows), div(sheetnumber - 1, nrows) + 1

function Base.getindex(p::BmPartition, i::Int)
    1 <= i <= p.ncols * p.nrows || throw(BoundsError(p, i))
    _SheetPartition(p, i)
end


function Base.getindex(bmp::BmPartition, r::Int, c::Int)
    1 <= r <= bmp.nrows || throw(BoundsError(bmp, (r, c)))
    1 <= c <= bmp.ncols || throw(BoundsError(bmp, (r, c)))
    i = r + (c - 1) * bmp.nrows
    _SheetPartition(bmp, i)
end
