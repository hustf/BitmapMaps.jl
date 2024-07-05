"""
    SheetBuilder

A SheetBuilder is one sheet in a collection defined by `SheetMatrixBuilder`. Direct construction is discouraged.

# Fields
- `pixel_origin_ref_to_bitmapmap`::Tuple{Int64, Int64}

The sheet's upper left corner given in a pixel coordinate system for the entire matrix of sheets. Intended for patching together sheets.

- `cell_iter`::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}

Iterate over every pixel / cell index `I` in the current sheet, starting top left. This field is the same for every sheet in a SheetMatrixBuilder.

- `f_I_to_utm`::Function

f_I_to_utm(I) transforms a cell index I for this sheet to an UTM coordinate. About UTM coordinates, see `GeoArrays.jl`.

- `sheet_number`::Int

Sequential index `i` of this sheet within the sheets matrix. You can get the sheet's row and column by calling `row_col_of_sheet(smb, i)`

- `density_pt_m⁻¹`::Int

When printing, 'DPI' means 'Dots Per Inch', which also is a density, but with inches instead of meters. Dots, cells, points and pixels is kind of the same thing.

- `pthsh`

The path to both input and output files for this sheet.

# Example

```
julia> using BitmapMaps

julia> smb = run_bitmapmap_pipeline(); # No matter how early the pipeline exits, 'smb' is generated.
...
...

julia> smb[2, 2]  # Let's examine one of the sheets
SheetBuilder(;pixel_origin_ref_to_bitmapmap = (2255, 3248),
                         cell_iter = CartesianIndices((1:3248, 1:2255)),
                        f_I_to_utm = (11638, 6928536),
                      sheet_number = 5,
                    density_pt_m⁻¹ = 11811,
                             pthsh = "bitmapmaps/default\\2 2  11638 6918792  18403 6928536")

julia> show_derived_properties(smb[2, 2]) # Let's see more human-readable details.

        [easting, northing] derived properties:
          Bounding Box (BB) SE-NW            = (11638 6918792)-(18403 6928536)
          Northeast internal corner          = (18400, 6928536) - most northeastern sample point
          Geo centre                         = (15020.5, 6.923664e6)
          Grid centre single                 = (15020, 6923665)
        Derived properties:
          Geographical area [km²]            = 66
          Adj. paper (width, height) [mm]    = (190.9, 275.0)
          Map scale                          = 1 : 35433 = 1 : (cell_to_utm_factor * density_pt_m⁻¹)
        External boundary box as Well Known Text (paste in e.g. https://nvdb-vegdata.github.io/nvdb-visrute/STM ):
          POLYGON ((11638 6918792, 18403 6918792, 18403 6928536, 11638 6928536, 11638 6918792))
```
"""
@kwdef struct SheetBuilder
    pixel_origin_ref_to_bitmapmap::Tuple{Int64, Int64}
    cell_iter::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}
    f_I_to_utm::Function
    sheet_number::Int
    density_pt_m⁻¹::Int
    pthsh::String
end

function SheetBuilder(pixel_origin_ref_to_bitmapmap, cell_iter, sheet_lower_left_utm, cell_to_utm_factor, sheet_number, density_pt_m⁻¹, pthsh)
    ny = size(cell_iter)[1]
    function f_I_to_utm(I::CartesianIndex)
        easting_offset =   (I[2] - 1) * cell_to_utm_factor
        northing_offset =  (ny - I[1] + 1) * cell_to_utm_factor
        (easting_offset, northing_offset) .+ sheet_lower_left_utm
    end
    SheetBuilder(;pixel_origin_ref_to_bitmapmap, cell_iter, f_I_to_utm, sheet_number, density_pt_m⁻¹, pthsh)
end

function Base.show(io::IO, ::MIME"text/plain", sb::SheetBuilder)
    print(io, "SheetBuilder(;")
    for sy in fieldnames(SheetBuilder)
        va = getfield(sb, sy)
        colwi = sy == :pixel_origin_ref_to_bitmapmap ? 9 : 34
        print(io, lpad("$sy", colwi), " = ")
        if sy == :f_I_to_utm
            print(io, sb.f_I_to_utm(first(sb.cell_iter)))
        else
            print(io, repr(va))
        end
        if sy !== :pthsh
            println(io, ", ")
        else
            print(io, ")")
        end
    end
end
function Base.show(io::IO, sb::SheetBuilder)
    print(io, "SheetBuilder(")
    print(io, sb.pixel_origin_ref_to_bitmapmap, ", ")
    print(io, sb.cell_iter.indices, ", ")
    print(io, "f(I) -> utm", ", ")
    print(io, sb.sheet_number)
    print(io, sb.density_pt_m⁻¹, ", ")
    print(io, repr(sb.pthsh), ", ")
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

A builder specifying a map divided in printable sheets. The sheet matrix builder is indexable and iterable as sheet builders.

Although direct construction is possible (see the test files), the intended construction is through calling `run_bitmapmap_pipeline`.

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
        colwi = sy == :southwest_corner ? 9 : 34
        print(io, lpad(repr(va), colwi))
        if sy !== :pth
            println(io, ", # ", "$sy")
        else
            print(io, ") # ", "$sy")
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
    opy = height_cell(smb) - sheet_height_cell(smb) - (r - 1) * sheet_height_cell(smb)
    pixel_origin_ref_to_bitmapmap = (opx, opy)
    sheet_lower_left_utm = smb.southwest_corner .+ smb.cell_to_utm_factor .* (opx, (r - 1) * sheet_height_cell(smb))
    # The local folder name includes row, column, min_easting, min_northing, max_easting, max_northing. This simplifies
    # the manual data ordering process (the web api requires user rights).
    min_easting, min_northing = sheet_lower_left_utm
    max_easting_external = min_easting + sheet_width_cell(smb) * smb.cell_to_utm_factor
    max_northing_external = min_northing + sheet_height_cell(smb) * smb.cell_to_utm_factor
    pthsh = joinpath(smb.pth, "$r $c  $min_easting $min_northing  $max_easting_external $max_northing_external")
    SheetBuilder(pixel_origin_ref_to_bitmapmap, cell_iter, sheet_lower_left_utm, smb.cell_to_utm_factor, sheet_number, smb.density_pt_m⁻¹, pthsh)
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
