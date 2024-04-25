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

- pix_iter::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}

Iterate over every pixel index `I` in the current sheet, starting top left. This field is the same for every sheet in a SheetMatrixBuilder.

- f_I_to_utm::Function

f_I_to_utm(I) transforms a pixel index I for this sheet to an UTM coordinate. About UTM coordinates, see `GeoArrays.jl`.

- sheet_number::Int

Sequential index `i` of this sheet. You can get the sheet's row and column by calling `row_col_of_sheet(smb, i)`


# Example

Note, this is not the recommended usage. Intended is: `run_bitmapmap_pipeline`.

```
julia> using BitmapMaps

julia> bitmapmap = let
           bm_pix_width, bm_pix_height = 9, 8              # Actual full bitmapmap sizes would be thousands of pixels
           sheet_pix_width, sheet_pix_height = 3, 4        # 12 sheets, over which bm_pixels can be evenly divided
           southwest_corner = (43999, 6909048)             # Utm zone coordinates, the corner is lower left in smb[1, 1] (defined below)
           pix_to_utm_factor = 2                           # One horizontal or vertical pixel distance equals 2 metres easting or 2 metres northing
           pth = "BitmapMaps/example"
           #
           SheetMatrixBuilder(bm_pix_width, bm_pix_height, sheet_pix_width, sheet_pix_height, southwest_corner, pix_to_utm_factor, pth)
       end;


julia> display.(bitmapmap);   # Broadcasting, mapping, iterations are defined.
       SheetBuilder(
               pixel_origin_ref_to_bitmapmap = (0, 4)
               pix_iter or I =                 CartesianIndices((1:4, 1:3))
               f_I_to_utm(first(pix_iter)) =   (43999, 6909054)
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
    map(sheetpartition.f_I_to_utm, sheetpartition.pix_iter)
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
    pix_iter::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}
    f_I_to_utm::Function
    sheet_number::Int
    pthsh::String
end

function SheetBuilder(sheet_lower_left_utm, pix_dist, pixel_origin_ref_to_bitmapmap, pix_iter, pthsh, sheet_number)
    sheet_pix_height = size(pix_iter)[1]
    function f_I_to_utm(I::CartesianIndex)
        easting_offset =   (I[2] - 1) * pix_dist
        northing_offset =  (sheet_pix_height - I[1] + 1) * pix_dist
        (easting_offset, northing_offset) .+ sheet_lower_left_utm
    end
    SheetBuilder(;pixel_origin_ref_to_bitmapmap, pix_iter, f_I_to_utm, pthsh, sheet_number)
end

function Base.show(io::IO, ::MIME"text/plain", shi::SheetBuilder)
    colwi = 32
    println(io, "SheetBuilder(")
    println(io, rpad("\tpixel_origin_ref_to_bitmapmap = ", colwi), shi.pixel_origin_ref_to_bitmapmap)
    println(io, rpad("\tpix_iter or I = ", colwi), shi.pix_iter)
    println(io, rpad("\tf_I_to_utm(first(pix_iter)) = ", colwi), shi.f_I_to_utm(first(shi.pix_iter)))
    println(io, rpad("\tpthsh = ", colwi), shi.pthsh)
    println(io, rpad("\tsheet_number = ", colwi), shi.sheet_number)
    println(io, " )")
end
function Base.show(io::IO, shi::SheetBuilder)
    print(io, "SheetBuilder(")
    print(io, shi.pixel_origin_ref_to_bitmapmap, ", ")
    print(io, shi.pix_iter.indices, ", ")
    print(io, "f(I) -> utm", ", ")
    print(io, "pthsh = ", shi.pthsh)
    print(io, shi.sheet_number)
    println(io, " )")
end



"""
    SheetMatrixBuilder(bm_pix_width, bm_pix_height, sheet_pix_width, sheet_pix_height, southwest_corner, pix_to_utm_factor, pth)
    SheetMatrixBuilder(pwi::Int, phe::Int, pdens_dpi::Int, nrc::Tuple{Int, Int}, southwest_corner::Tuple{Int, Int}, pix_to_utm_factor::Int, pth)

An object specifying the map and helping with making it.

See `pipeline.jl/define_job(;kwds)` for typical construction.

See `SheetBuilder` and `resource/matrix_sheet_pix_utm.svg` for an example.
"""
struct SheetMatrixBuilder
    bm_pix_width::Int
    bm_pix_height::Int
    sheet_pix_width::Int
    sheet_pix_height::Int
    nrows::Int
    ncols::Int
    southwest_corner::Tuple{Int, Int}
    sheet_indices::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}
    pix_to_utm_factor::Int
    pth::String
    function SheetMatrixBuilder(bm_pix_width, bm_pix_height, sheet_pix_width, sheet_pix_height, southwest_corner, pix_to_utm_factor, pth)
        ncols = bm_pix_width / sheet_pix_width |> ceil |> Integer
        nrows = bm_pix_height / sheet_pix_height |> ceil |> Integer
        if ncols < 1 || nrows < 1
            throw(ArgumentError("bm_pix_... must be greater than corresponding sheet_pix... "))
        end
        if sheet_pix_width * ncols !== bm_pix_width 
            throw(ArgumentError("bm_pix_width $bm_pix_width must be an integer times sheet_pix_width $sheet_pix_width"))
        end
        if sheet_pix_height * nrows !== bm_pix_height 
            throw(ArgumentError("bm_pix_height $bm_pix_height must be an integer times sheet_pix_height $sheet_pix_width"))
        end
        sheet_indices = CartesianIndices((1:nrows, 1:ncols))
        new(bm_pix_width, bm_pix_height, sheet_pix_width, sheet_pix_height, nrows, ncols, southwest_corner, sheet_indices, pix_to_utm_factor, pth)
    end
end

function SheetMatrixBuilder(pwi::Int, phe::Int, pdens_dpi::Int, nrc::Tuple{Int, Int}, southwest_corner::Tuple{Int, Int}, pix_to_utm_factor::Int, pth::String)
    sheet_pix_width = floor(pwi * pdens_dpi / 25.4)
    sheet_pix_height = floor(phe * pdens_dpi / 25.4)
    nrows, ncols = nrc
    bm_pix_width = ncols * sheet_pix_width
    bm_pix_height = nrows * sheet_pix_height
    SheetMatrixBuilder(bm_pix_width, bm_pix_height, sheet_pix_width, sheet_pix_height, southwest_corner, pix_to_utm_factor, pth)
end

function Base.show(io::IO, ::MIME"text/plain", p::SheetMatrixBuilder)
    colwi = 32
    println(io, "SheetMatrixBuilder(")
    for sy in fieldnames(SheetMatrixBuilder)
        va = getfield(p, sy)
        print(io, rpad("\t$sy", colwi), " = ", va)
        println(sy !== :pix_to_utm_factor ? "" : ")")
    end
end


function _SheetPartition(p::SheetMatrixBuilder, sheet_number::Int)
    # No bounds in this "private" function, 'iterate' and 'getindex' should do that.
    pix_iter = CartesianIndices((1:p.sheet_pix_height, 1:p.sheet_pix_width))
    r, c = row_col_of_sheet(p, sheet_number)
    opx = (c - 1) * p.sheet_pix_width
    opy = p.bm_pix_height - p.sheet_pix_height - (r - 1) * p.sheet_pix_height
    pixel_origin_ref_to_bitmapmap = (opx, opy)
    sheet_lower_left_utm = p.southwest_corner .+ p.pix_to_utm_factor .* (opx, (r - 1) * p.sheet_pix_height)
    # The local folder name includes row, column, min_easting, min_northing, max_easting, max_northing. This simplifies
    # the manual data ordering process (the web api requires user rights).
    min_easting, min_northing = sheet_lower_left_utm
    max_easting_external = min_easting + p.sheet_pix_width * p.pix_to_utm_factor
    max_northing_external = min_northing + p.sheet_pix_height * p.pix_to_utm_factor
    pthsh = joinpath(p.pth, "$r $c  $min_easting $min_northing  $max_easting_external $max_northing_external")
    SheetBuilder(sheet_lower_left_utm, p.pix_to_utm_factor, pixel_origin_ref_to_bitmapmap, pix_iter, pthsh, sheet_number)
end

Base.iterate(smb::SheetMatrixBuilder) = _SheetPartition(smb, 1), _SheetPartition(smb, 2)

function Base.iterate(p::SheetMatrixBuilder, state::SheetBuilder)
    state.sheet_number > p.nrows * p.ncols && return nothing
    # return: the one to use next, the one after that
    state, _SheetPartition(p, state.sheet_number + 1)
end


Base.length(smb::SheetMatrixBuilder) = smb.nrows * smb.ncols
Base.size(smb::SheetMatrixBuilder) = (smb.nrows, smb.ncols)
Base.axes(smb::SheetMatrixBuilder, d) = axes(smb)[d]
Base.axes(smb::SheetBuilder) = map(Base.oneto, size(smb))



row_col_of_sheet(p::SheetMatrixBuilder, sheet_number::Int) = row_col_of_sheet(p.nrows, sheet_number)
row_col_of_sheet(nrows::Int, sheetnumber::Int) = mod1(sheetnumber, nrows), div(sheetnumber - 1, nrows) + 1

function Base.getindex(p::SheetMatrixBuilder, i::Int)
    1 <= i <= p.ncols * p.nrows || throw(BoundsError(p, i))
    _SheetPartition(p, i)
end


function Base.getindex(smb::SheetMatrixBuilder, r::Int, c::Int)
    1 <= r <= smb.nrows || throw(BoundsError(smb, (r, c)))
    1 <= c <= smb.ncols || throw(BoundsError(smb, (r, c)))
    i = r + (c - 1) * smb.nrows
    _SheetPartition(smb, i)
end

#=
# Utilty functions

northeast_internal_corner(p::SheetBuilder) = p.f_I_to_utm(p.pix_iter[1, end])
northeast_internal_corner(p::SheetMatrixBuilder) = p.southwest_corner .+ p.pix_to_utm_factor .* (p.bm_pix_width - 1, p.bm_pix_height - 1)
function northeast_external_corner(p::SheetBuilder)
    pixel_distance_meter = p.f_I_to_utm(p.pix_iter[1])[2] - p.f_I_to_utm(p.pix_iter[2])[2]
    northeast_internal_corner(p) .+ (pixel_distance_meter, pixel_distance_meter)
end
northeast_external_corner(p::SheetMatrixBuilder) = p.southwest_corner .+ p.pix_to_utm_factor .* (p.bm_pix_width, p.bm_pix_height)
southwest_external_corner(p::SheetBuilder) = p.f_I_to_utm(p.pix_iter[end, 1]) .+ (p.f_I_to_utm(p.pix_iter[end, 1]) .- p.f_I_to_utm(p.pix_iter[end - 1, 1]))
southwest_external_corner(p::SheetMatrixBuilder) = p.southwest_corner
southwest_internal_corner(p::SheetBuilder) = p.f_I_to_utm(p.pix_iter[end, 1]) .+ (p.f_I_to_utm(p.pix_iter[end, 1]) .- p.f_I_to_utm(p.pix_iter[end - 1, 1]))

northwest_corner()

"""
    southwest_corner(p::SheetBuilder)
    southwest_corner(p::SheetMatrixBuilder)
    ---> Tuple{Int, Int}

The utm position is included in this partition. The geographical sample point lies one pixel distance north of the corner.
"""

"""
    northeast_external_corner(p::SheetBuilder)
    northeast_external_corner(p::SheetMatrixBuilder)
    ---> Tuple{Int, Int}

This utm position is NOT part of the sheet or bitmap, but equals the southwest position of the next one. See 'northeast_external_corner'.
"""

"""
    northeast_internal_corner(p::SheetBuilder)
    northeast_internal_corner(p::SheetMatrixBuilder)
    ---> Tuple{Int, Int}

This utm position is part of the sheet or bitmap. See 'northeast_external_corner'
"""







"""
    geo_grid_centre_single(p)
    ---> Tuple{Int, Int}

When the number of sample points is even, the actual centre lies between elements.
This returns the most easterly and southerly elements in those cases.

By 'the grid' we mean internal sample points of p.

Also see 'geometric_grid_centre'.
"""
geo_grid_centre_single(p) = div.(southwest_corner(p) .+ northeast_internal_corner(p), 2)


"""
    geo_centre(p)
    ---> Tuple{Float64, Float64}

This returns the centre point as a continuous coordinate. Since corner coordinates of
`p` are integers, the centre lies at 'n.0' for odd dimensions of `p`, or 'n.5' for even.

Also see 'geo_grid_centre_single'.
"""
geo_centre(p) = (southwest_corner(p) .+ northeast_external_corner(p)) ./ 2

bounding_box_external_string(p) = replace("$(southwest_corner(p))-$(northeast_external_corner(p))", "," => "")



function geo_area(p)
    s, w = southwest_corner(p)
    n, e = northeast_external_corner(p)
    (n - s) * (e - w)
end

# Well known text (for pasting elsewhere)

function bounding_box_closed_polygon_string(p::SheetBuilder)
    x1, y1 = southwest_corner(p)
    x3, y3 = northeast_external_corner(p)
    x2, y2 = x3, y1
    x4, y4 = x1, y3
   "($x1 $y1, $x2 $y2, $x3 $y3, $x4 $y4, $x1 $y1)"
end

function bounding_box_closed_polygon_string(p::SheetMatrixBuilder)
    s = ""
    for shp in p
        s *= "$(bounding_box_closed_polygon_string(shp))"
        if shp.sheet_number !== length(p)
            s *= ",\n\t\t"
        end
    end
    s
end

bounding_box_polygon_string(p) = "POLYGON ($(bounding_box_closed_polygon_string(p)))"

function show_augmented(p::SheetMatrixBuilder)
    printstyled("Bitmapmap configuration based on .ini file  ", color = :green, bold=:true)
    show(stdout, MIME("text/plain"), p)
    show_augmented_properties(p)
end

function show_augmented_properties(p)
    printstyled("\tAugmented properties (all as (easting, northing)): \n", color = :green)
    println("\t  ", rpad("Geo centre = ",        35), geo_centre(p))
    println("\t  ", rpad("Grid centre single = ", 35), geo_grid_centre_single(p))
    println("\t  ", rpad("Northeast external corner = ",     35), northeast_external_corner(p))
    println("\t  ", rpad("Northeast internal corner = ",     35), northeast_internal_corner(p), " - most northeastern sample point")
    println("\t  ", rpad("Bounding Box (BB) SE-NW = ", 35), bounding_box_external_string(p))
    if p isa SheetMatrixBuilder
        println("\t  ", rpad("Geographical area [km²] = ", 35), Int(round(geo_area(p) / 1e6)), "         Per sheet: ", round(geo_area(first(p)) / 1e6, digits = 1), "  km²   Single file export limit: 16 km²")
    else
        println("\t  ", rpad("Geographical area [km²] = ", 35), Int(round(geo_area(p) / 1e6)))
    end
    printstyled("\tBBs of sheets as Well Known Text ", color = :green)
    printstyled("(paste in e.g. https://nvdb-vegdata.github.io/nvdb-visrute/STM ):\n", color = :light_black)
    println("\t  ", bounding_box_polygon_string(p))
end
=#