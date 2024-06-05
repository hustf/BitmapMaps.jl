#=
Inspired by Luxor.Partition, adapted so as to ensure that
 - an integer number of sheets covers the entire area
 - the oposite direction of 'northing' and 'image matrix y+' does not cause inaccuracies (which it might do).
 - Hide away the bug nest of sampling distances.
=#

"""
    SheetBuilder

A SheetBuilder is one sheet in a collection defined by `SheetMatrixBuilder`.

# Fields
- pixel_origin_ref_to_bitmapmap::Tuple{Int64, Int64}

The sheet's upper left corner given in a pixel coordinate system for the entire sheets matrix.

- cell_iter::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}

Iterate over every pixel / cell index `I` in the current sheet, starting top left. This field is the same for every sheet in a SheetMatrixBuilder.

- f_I_to_utm::Function

f_I_to_utm(I) transforms a cell index I for this sheet to an UTM coordinate. About UTM coordinates, see `GeoArrays.jl`.

- sheet_number::Int

Sequential index `i` of this sheet. You can get the sheet's row and column by calling `row_col_of_sheet(smb, i)`


# Example

Note, this is not the recommended usage. Intended is: `run_bitmapmap_pipeline`.

```
julia> using BitmapMaps

julia> smb = let
           bm_cell_width, bm_cell_height = 9, 8              # Actual full bitmapmap sizes would be thousands of pixels
           sheet_cell_width, sheet_cell_height = 3, 4        # 12 sheets, over which bm_pixels can be evenly divided
           southwest_corner = (43999, 6909048)             # Utm zone coordinates, the corner is lower left in smb[1, 1] (defined below)
           cell_to_utm_factor = 2                           # One horizontal or vertical cell distance equals 2 metres easting or 2 metres northing
           sheet_cell_width_mm = 191                        # This value is the default in the .ini file.
           pth = "BitmapMaps/example"
           #
           SheetMatrixBuilder(bm_cell_width, bm_cell_height, sheet_cell_width, sheet_cell_height, southwest_corner, cell_to_utm_factor, sheet_cell_width_mm, pth)
       end;

TODO: UPDATE
julia> display.(smg;   # Broadcasting, mapping, iterations are defined.
       SheetBuilder(
               pixel_origin_ref_to_bitmapmap = (0, 4)
               cell_iter or I =                 CartesianIndices((1:4, 1:3))
               f_I_to_utm(first(cell_iter)) =   (43999, 6909054)
               sheet_number =                  1
        )
        ⋮
               sheet_number =                  12
        )

julia> bitmapmap[6] |> println     # Single index, short form printing.
SheetBuilder((6, 0), (1:4, 1:3), f(I) -> utm, 6 )

julia> bitmapmap[2, 3] |> println  # Row, col index, same object
SheetBuilder((6, 0), (1:4, 1:3), f(I) -> utm, 6 )

julia> utm_coordinates_of_pixels_in_sheet = let
    sheetpartition = bitmapmap[1]
    map(sheetpartition.f_I_to_utm, sheetpartition.cell_iter)
end
4×3 Matrix{Tuple{Int64, Int64}}:
(43999, 6909054)  (44001, 6909054)  (44003, 6909054)
(43999, 6909052)  (44001, 6909052)  (44003, 6909052)
(43999, 6909050)  (44001, 6909050)  (44003, 6909050)
(43999, 6909048)  (44001, 6909048)  (44003, 6909048)
```
"""
@kwdef struct SheetBuilder
    pixel_origin_ref_to_bitmapmap::Tuple{Int64, Int64}
    cell_iter::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}
    f_I_to_utm::Function
    sheet_number::Int
    sheet_width_mm::Int
    pthsh::String
end

function SheetBuilder(sheet_lower_left_utm, cell_size, pixel_origin_ref_to_bitmapmap, cell_iter, sheet_width_mm, pthsh, sheet_number)
    sheet_cell_height = size(cell_iter)[1]
    function f_I_to_utm(I::CartesianIndex)
        easting_offset =   (I[2] - 1) * cell_size
        northing_offset =  (sheet_cell_height - I[1] + 1) * cell_size
        (easting_offset, northing_offset) .+ sheet_lower_left_utm
    end
    SheetBuilder(;pixel_origin_ref_to_bitmapmap, cell_iter, f_I_to_utm, sheet_width_mm, pthsh, sheet_number)
end

function Base.show(io::IO, ::MIME"text/plain", sb::SheetBuilder)
    colwi = 32
    println(io, "SheetBuilder(")
    println(io, rpad("    pixel_origin_ref_to_bitmapmap = ", colwi), sb.pixel_origin_ref_to_bitmapmap)
    println(io, rpad("    cell_iter or I = ", colwi), sb.cell_iter)
    println(io, rpad("    f_I_to_utm(first(cell_iter)) = ", colwi), sb.f_I_to_utm(first(sb.cell_iter)))
    println(io, rpad("    sheet_width_mm = ", colwi), sb.sheet_width_mm)
    println(io, rpad("    pthsh = ", colwi), sb.pthsh)
    println(io, rpad("    sheet_number = ", colwi), sb.sheet_number)
    println(io, " )")
end
function Base.show(io::IO, sb::SheetBuilder)
    print(io, "SheetBuilder(")
    print(io, sb.pixel_origin_ref_to_bitmapmap, ", ")
    print(io, sb.cell_iter.indices, ", ")
    print(io, "f(I) -> utm", ", ")
    print(io, sb.sheet_width_mm, ", ")
    print(io, repr(sb.pthsh), ", ")
    print(io, sb.sheet_number)
    println(io, " )")
end



"""
    SheetMatrixBuilder(bm_cell_width, bm_cell_height, sheet_cell_width, sheet_cell_height, southwest_corner, cell_to_utm_factor, pth)
    SheetMatrixBuilder(pwi::Int, phe::Int, pdens_dpi::Int, nrc::Tuple{Int, Int}, southwest_corner::Tuple{Int, Int}, cell_to_utm_factor::Int, pth)
    SheetMatrixBuilder(bm_cell_width::Int, bm_cell_height::Int, sheet_cell_width::Int, sheet_cell_height::Int, nrows::Int, ncols::Int, southwest_corner::Tuple{Int, Int}, 
        sheet_indices::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}, cell_to_utm_factor::Int, pth::String)

An object specifying the map and helping with making it.
These constructors aren't intended for use directly, but should allow parseable output. Instead, see `pipeline.jl/define_builder(;kwds)` for typical construction.

See `SheetBuilder` and `resource/matrix_sheet_cell_utm.svg` for an example.
"""
struct SheetMatrixBuilder
    bm_cell_width::Int
    bm_cell_height::Int
    sheet_cell_width::Int
    sheet_cell_height::Int
    nrows::Int
    ncols::Int
    southwest_corner::Tuple{Int, Int}
    sheet_indices::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}
    cell_to_utm_factor::Int
    sheet_width_mm::Int
    pth::String
    function SheetMatrixBuilder(bm_cell_width, bm_cell_height, sheet_cell_width, sheet_cell_height, southwest_corner, cell_to_utm_factor, sheet_width_mm, pth)
        ncols = bm_cell_width / sheet_cell_width |> ceil |> Integer
        nrows = bm_cell_height / sheet_cell_height |> ceil |> Integer
        if ncols < 1 || nrows < 1
            throw(ArgumentError("bm_cell_... must be greater than corresponding sheet_pix... "))
        end
        if sheet_cell_width * ncols !== bm_cell_width
            throw(ArgumentError("bm_cell_width $bm_cell_width must be an integer times sheet_cell_width $sheet_cell_width"))
        end
        if sheet_cell_height * nrows !== bm_cell_height
            throw(ArgumentError("bm_cell_height $bm_cell_height must be an integer times sheet_cell_height $sheet_cell_width"))
        end
        sheet_indices = CartesianIndices((1:nrows, 1:ncols))
        new(bm_cell_width, bm_cell_height, sheet_cell_width, sheet_cell_height, nrows, ncols, southwest_corner, sheet_indices, cell_to_utm_factor, sheet_width_mm, pth)
    end
end

function SheetMatrixBuilder(pwi::Int, phe::Int, pdens_dpi::Int, nrc::Tuple{Int, Int}, southwest_corner::Tuple{Int, Int}, cell_to_utm_factor::Int, pth::String)
    sheet_cell_width = floor(pwi * pdens_dpi / 25.4)
    sheet_cell_height = floor(phe * pdens_dpi / 25.4)
    nrows, ncols = nrc
    bm_cell_width = ncols * sheet_cell_width
    bm_cell_height = nrows * sheet_cell_height
    sheet_width_mm = Int(ceil(sheet_cell_width / (pdens_dpi / 25.4)))
    SheetMatrixBuilder(bm_cell_width, bm_cell_height, sheet_cell_width, sheet_cell_height, southwest_corner, cell_to_utm_factor, sheet_width_mm, pth)
end
# This method is 'over-determined', but allows parsing output (that also contains the same information more than twice)
function SheetMatrixBuilder(bm_cell_width::Int, bm_cell_height::Int, sheet_cell_width::Int, sheet_cell_height::Int, nrows::Int, ncols::Int, southwest_corner::Tuple{Int, Int}, 
    sheet_indices::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}, cell_to_utm_factor::Int, sheet_width_mm::Int, pth::String)
    @assert nrows == Int(ceil(bm_cell_height / sheet_cell_height))
    @assert ncols == Int(ceil(bm_cell_width / sheet_cell_width ))
    @assert sheet_indices == CartesianIndices((1:nrows, 1:ncols))
    SheetMatrixBuilder(bm_cell_width, bm_cell_height, sheet_cell_width, sheet_cell_height, southwest_corner, cell_to_utm_factor, sheet_width_mm, pth)
end



function Base.show(io::IO, ::MIME"text/plain", smb::SheetMatrixBuilder)
    print(io, "SheetMatrixBuilder(")
    for sy in fieldnames(SheetMatrixBuilder)
        va = getfield(smb, sy)
        colwi = sy == :bm_cell_width ? 9 : 28
        print(io, lpad(repr(va), colwi))
        if sy !== :pth
            println(io, ", # ", "$sy")
        else
            println(io, ") # ", "$sy")
        end
    end
end


function _SheetBuilder(smb::SheetMatrixBuilder, sheet_number::Int)
    # No bounds in this "private" function, 'iterate' and 'getindex' should do that.
    cell_iter = CartesianIndices((1:smb.sheet_cell_height, 1:smb.sheet_cell_width))
    r, c = row_col_of_sheet(smb, sheet_number)
    opx = (c - 1) * smb.sheet_cell_width
    opy = smb.bm_cell_height - smb.sheet_cell_height - (r - 1) * smb.sheet_cell_height
    pixel_origin_ref_to_bitmapmap = (opx, opy)
    sheet_lower_left_utm = smb.southwest_corner .+ smb.cell_to_utm_factor .* (opx, (r - 1) * smb.sheet_cell_height)
    sheet_width_mm = smb.sheet_width_mm
    # The local folder name includes row, column, min_easting, min_northing, max_easting, max_northing. This simplifies
    # the manual data ordering process (the web api requires user rights).
    min_easting, min_northing = sheet_lower_left_utm
    max_easting_external = min_easting + smb.sheet_cell_width * smb.cell_to_utm_factor
    max_northing_external = min_northing + smb.sheet_cell_height * smb.cell_to_utm_factor
    pthsh = joinpath(smb.pth, "$r $c  $min_easting $min_northing  $max_easting_external $max_northing_external")
    SheetBuilder(sheet_lower_left_utm, smb.cell_to_utm_factor, pixel_origin_ref_to_bitmapmap, cell_iter, sheet_width_mm, pthsh, sheet_number)
end

Base.iterate(smb::SheetMatrixBuilder) = _SheetBuilder(smb, 1), _SheetBuilder(smb, 2)

function Base.iterate(smb::SheetMatrixBuilder, state::SheetBuilder)
    state.sheet_number > smb.nrows * smb.ncols && return nothing
    # return: the one to use next, the one after that
    state, _SheetBuilder(smb, state.sheet_number + 1)
end


Base.length(smb::SheetMatrixBuilder) = smb.nrows * smb.ncols
Base.size(smb::SheetMatrixBuilder) = (smb.nrows, smb.ncols)
Base.axes(smb::SheetMatrixBuilder, d) = axes(smb)[d]
Base.axes(smb::SheetBuilder) = map(Base.oneto, size(smb))



row_col_of_sheet(smb::SheetMatrixBuilder, sheet_number::Int) = row_col_of_sheet(smb.nrows, sheet_number)
row_col_of_sheet(nrows::Int, sheetnumber::Int) = mod1(sheetnumber, nrows), div(sheetnumber - 1, nrows) + 1

function Base.getindex(smb::SheetMatrixBuilder, i::Int)
    1 <= i <= smb.ncols * smb.nrows || throw(BoundsError(smb, i))
    _SheetBuilder(smb, i)
end


function Base.getindex(smb::SheetMatrixBuilder, r::Int, c::Int)
    1 <= r <= smb.nrows || throw(BoundsError(smb, (r, c)))
    1 <= c <= smb.ncols || throw(BoundsError(smb, (r, c)))
    i = r + (c - 1) * smb.nrows
    _SheetBuilder(smb, i)
end
