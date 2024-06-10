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
           bm_width_cell, bm_height_cell = 9, 8              # Actual full bitmapmap sizes would be thousands of pixels
           sheet_width_cell, sheet_height_cell = 3, 4        # 12 sheets, over which bm_pixels can be evenly divided
           southwest_corner = (43999, 6909048)             # Utm zone coordinates, the corner is lower left in smb[1, 1] (defined below)
           cell_to_utm_factor = 2                           # One horizontal or vertical cell distance equals 2 metres easting or 2 metres northing
           sheet_cell_width_mm = 191                        # This value is the default in the .ini file.
           pth = "BitmapMaps/example"
           #
           SheetMatrixBuilder(bm_width_cell, bm_height_cell, sheet_width_cell, sheet_height_cell, southwest_corner, cell_to_utm_factor, sheet_cell_width_mm, pth)
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
    sheet_height_cell = size(cell_iter)[1]
    function f_I_to_utm(I::CartesianIndex)
        easting_offset =   (I[2] - 1) * cell_size
        northing_offset =  (sheet_height_cell - I[1] + 1) * cell_size
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
    SheetMatrixBuilder(southwest_corner::Tuple{Int, Int},
        sheet_indices::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}},
        cell_to_utm_factor::Int,
        sheet_width_mm::Int,
        sheet_height_mm::Int,
        density_pt_m⁻¹::Int,
        pth::String)
    SheetMatrixBuilder(southwest_corner::Tuple{Int, Int},
        nrc::Tuple{Int, Int},
        cell_to_utm_factor::Int,
        sheet_width_mm::Int,
        sheet_height_mm::Int,
        density_pt_m⁻¹::Int,
        pth::String)

A builder specifying a map divided in printable sheets. This builder is indexable as sheet builders.
These constructors aren't intended for use directly. Instead, see `pipeline.jl/define_builder(;kwds)` for typical construction.

See `SheetBuilder` and `resource/matrix_sheet_cell_utm.svg` for an example.
"""
struct SheetMatrixBuilder
    southwest_corner::Tuple{Int, Int}
    sheet_indices::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}
    cell_to_utm_factor::Int
    sheet_width_mm::Int
    sheet_height_mm::Int
    density_pt_m⁻¹::Int
    pth::String
    function SheetMatrixBuilder(southwest_corner, sheet_indices::CartesianIndices, cell_to_utm_factor, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, pth)
        ! isempty(sheet_indices) || throw(ArgumentError("sheet_indices"))
        cell_to_utm_factor > 0 || throw(ArgumentError("cell_to_utm_factor"))
        sheet_width_mm > 0 || throw(ArgumentError("sheet_width_mm "))
        sheet_height_mm > 0 || throw(ArgumentError("sheet_height_mm"))
        density_pt_m⁻¹ > 0 || throw(ArgumentError("density_pt_m⁻¹"))
        ! isempty(pth) || throw(ArgumentError("pth"))
        new(southwest_corner, sheet_indices, cell_to_utm_factor, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, pth)
    end
end
function SheetMatrixBuilder(southwest_corner::Tuple{Int, Int},
    nrc::Tuple{Int, Int},
    cell_to_utm_factor::Int,
    sheet_width_mm::Int,
    sheet_height_mm::Int,
    density_pt_m⁻¹::Int,
    pth::String)
    # Make the sheet indices
    nrows, ncols = nrc
    sheet_indices = CartesianIndices((1:nrows, 1:ncols))
    SheetMatrixBuilder(southwest_corner, sheet_indices, cell_to_utm_factor, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, pth)
end

function Base.show(io::IO, ::MIME"text/plain", smb::SheetMatrixBuilder)
    print(io, "SheetMatrixBuilder(")
    for sy in fieldnames(SheetMatrixBuilder)
        va = getfield(smb, sy)
        colwi = sy == :bm_width_cell ? 9 : 28
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
    nx = sheet_width_cell(smb)
    ny = sheet_height_cell(smb)
    cell_iter = CartesianIndices((1:ny, 1:nx))
    r, c = row_col_of_sheet(smb, sheet_number)
    opx = (c - 1) * sheet_width_cell(smb)
    opy = bm_height_cell(smb) - sheet_height_cell(smb) - (r - 1) * sheet_height_cell(smb)
    pixel_origin_ref_to_bitmapmap = (opx, opy)
    sheet_lower_left_utm = smb.southwest_corner .+ smb.cell_to_utm_factor .* (opx, (r - 1) * sheet_height_cell(smb))
    sheet_width_mm = smb.sheet_width_mm
    # The local folder name includes row, column, min_easting, min_northing, max_easting, max_northing. This simplifies
    # the manual data ordering process (the web api requires user rights).
    min_easting, min_northing = sheet_lower_left_utm
    max_easting_external = min_easting + sheet_width_cell(smb) * smb.cell_to_utm_factor
    max_northing_external = min_northing + sheet_height_cell(smb) * smb.cell_to_utm_factor
    pthsh = joinpath(smb.pth, "$r $c  $min_easting $min_northing  $max_easting_external $max_northing_external")
    SheetBuilder(sheet_lower_left_utm, smb.cell_to_utm_factor, pixel_origin_ref_to_bitmapmap, cell_iter, sheet_width_mm, pthsh, sheet_number)
end

Base.iterate(smb::SheetMatrixBuilder) = _SheetBuilder(smb, 1), _SheetBuilder(smb, 2)

function Base.iterate(smb::SheetMatrixBuilder, state::SheetBuilder)
    state.sheet_number > nrows(smb) * ncols(smb) && return nothing
    # return: the one to use next, the one after that
    state, _SheetBuilder(smb, state.sheet_number + 1)
end


Base.length(smb::SheetMatrixBuilder) = nrows(smb) * ncols(smb)
Base.size(smb::SheetMatrixBuilder) = (nrows(smb), ncols(smb))
Base.axes(smb::SheetMatrixBuilder, d) = axes(smb)[d]
Base.axes(smb::SheetBuilder) = map(Base.oneto, size(smb))



row_col_of_sheet(smb::SheetMatrixBuilder, sheet_number::Int) = row_col_of_sheet(nrows(smb), sheet_number)
row_col_of_sheet(nrows::Int, sheetnumber::Int) = mod1(sheetnumber, nrows), div(sheetnumber - 1, nrows) + 1

function Base.getindex(smb::SheetMatrixBuilder, i::Int)
    1 <= i <= ncols(smb) * nrows(smb) || throw(BoundsError(smb, i))
    _SheetBuilder(smb, i)
end


function Base.getindex(smb::SheetMatrixBuilder, r::Int, c::Int)
    1 <= r <= nrows(smb) || throw(BoundsError(smb, (r, c)))
    1 <= c <= ncols(smb) || throw(BoundsError(smb, (r, c)))
    i = r + (c - 1) * nrows(smb)
    _SheetBuilder(smb, i)
end
