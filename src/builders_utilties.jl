# Utilty functions
# Argument names
# - smb::SheetMatrixBuilder
# - sb::SheetBuilder
# - p - duck typed, can be a GeoArray or a SheetBuilder.
#
# A few methods are extended for GeoArray in `geoarray_utilties`


sheet_width_cell(smb::SheetMatrixBuilder) = Int(floor(smb.sheet_width_mm * smb.density_pt_m⁻¹ / 1000))
sheet_height_cell(smb::SheetMatrixBuilder) = Int(floor(smb.sheet_height_mm * smb.density_pt_m⁻¹ / 1000))

nrows(smb::SheetMatrixBuilder) = first(size(smb.sheet_indices))
nrows(sb::SheetBuilder) = first(size(sb.cell_iter))
ncols(smb::SheetMatrixBuilder) = last(size(smb.sheet_indices))
ncols(sb::SheetBuilder) = last(size(sb.cell_iter))

width_cell(smb::SheetMatrixBuilder) = ncols(smb) * sheet_width_cell(smb)
width_cell(sb::SheetBuilder) = ncols(sb)
height_cell(smb::SheetMatrixBuilder) = nrows(smb) * sheet_height_cell(smb)
height_cell(sb::SheetBuilder) = nrows(sb)

cell_to_utm_factor(sb::SheetBuilder) = sb.f_I_to_utm(CartesianIndex(1, 2))[1] - sb.f_I_to_utm(CartesianIndex(1, 1))[1]
map_scale_denominator(smb::SheetMatrixBuilder) = map_scale_denominator(first(smb))
map_scale_denominator(sb::SheetBuilder) = sb.density_pt_m⁻¹ * cell_to_utm_factor(sb)

"""
    width_adjusted_mm(smb::SheetMatrixBuilder)
    ---> Float64

Adjusted width per sheet is smaller than or equal to `sheet_width_mm` due to an integer number of cells.
"""
width_adjusted_mm(smb::SheetMatrixBuilder) = ncols(smb) * width_adjusted_mm(first(smb))
width_adjusted_mm(sb::SheetBuilder) = 1000 * width_cell(sb) / sb.density_pt_m⁻¹

"""
height_adjusted_mm(smb::SheetMatrixBuilder)
---> Float64

Adjusted height per sheet  is smaller than or equal to `sheet_height_mm` due to an integer number of cells.
"""
height_adjusted_mm(smb::SheetMatrixBuilder) = ncols(smb) * height_adjusted_mm(first(smb))
height_adjusted_mm(sb::SheetBuilder) = 1000 * height_cell(sb) / sb.density_pt_m⁻¹


"""
    southwest_corner(smb::SheetMatrixBuilder)
    northeast_corner(smb::SheetMatrixBuilder)
    northwest_corner(smb::SheetMatrixBuilder)
    southeast_corner(smb::SheetMatrixBuilder)
    northwest_corner(p::SheetBuilder)
    southeast_internal_corner(p::SheetBuilder)
    southeast_external_corner(p::SheetBuilder)
    southwest_internal_corner(p::SheetBuilder)
    southwest_external_corner(p::SheetBuilder)
    northeast_internal_corner(p::SheetBuilder)
    northeast_external_corner(p::SheetBuilder)
    ---> Tuple(Int, Int)

An entire SheetMatrixBuilder is divided into sheets. Every SheetMatrixBuilder method refers to 'external' borders.
A SheetBuilder is a subset of a SheetMatrixBuilder. It shares corners with neighbouring sheets.
Each sheet is divided into cells (i.e. pixels), and each cell is represented by a single geographical position or sampling point.

'external' refers to outside dimensions of the sheet or BitMap, which are larger than 'internal'.
'internal' refers to the last sampling point within or on the border of 'external' borders.

Each geographical sampling point corresponds to the top left of its cell. Hence, 'external' and 'internal'
are indentical for the northwest corner of a sheet.
"""
southwest_corner(smb::SheetMatrixBuilder) = smb.southwest_corner
"ref. southwest_corner"
northeast_corner(smb::SheetMatrixBuilder) = southwest_corner(smb) .+ smb.cell_to_utm_factor .* (width_cell(smb), height_cell(smb))
"ref. southwest_corner"
northwest_corner(smb::SheetMatrixBuilder) = southwest_corner(smb) .+ smb.cell_to_utm_factor .* (0, height_cell(smb))
"ref. southwest_corner"
southeast_corner(smb::SheetMatrixBuilder) = southwest_corner(smb) .+ smb.cell_to_utm_factor .* (width_cell(smb), 0)
"ref. southwest_corner"
northwest_corner(sb::SheetBuilder) = sb.f_I_to_utm(CartesianIndex(1, 1))
"ref. southwest_corner"
southwest_external_corner(smb::SheetMatrixBuilder) = southwest_corner(smb)
"ref. southwest_corner"
northeast_external_corner(smb::SheetMatrixBuilder) = northeast_corner(smb)
"ref. southwest_corner"
northwest_external_corner(smb::SheetMatrixBuilder) = northwest_corner(smb)
"ref. southwest_corner"
southeast_external_corner(smb::SheetMatrixBuilder) = southeast_corner(smb)
"ref. southwest_corner"
southeast_internal_corner(p::SheetBuilder) = p.f_I_to_utm(p.cell_iter[end, end])
"ref. southwest_corner"
southeast_internal_corner(smb::SheetMatrixBuilder) = southeast_internal_corner(smb[1, end])
"ref. southwest_corner"
southeast_external_corner(p::SheetBuilder) = 2 .* southeast_internal_corner(p) .- p.f_I_to_utm(p.cell_iter[end - 1, end - 1])
"ref. southwest_corner"
southwest_internal_corner(p) = (northwest_corner(p)[1], southeast_internal_corner(p)[2])
southwest_internal_corner(smb::SheetMatrixBuilder) = southwest_internal_corner(smb[1, 1])
"ref. southwest_corner"
southwest_external_corner(p) = (northwest_corner(p)[1], southeast_external_corner(p)[2])
"ref. southwest_corner"
northeast_internal_corner(p) = (southeast_internal_corner(p)[1], northwest_corner(p)[2])
northeast_internal_corner(smb::SheetMatrixBuilder) = northeast_internal_corner(smb[end, end])
"ref. southwest_corner"
northeast_external_corner(p) = (southeast_external_corner(p)[1], northwest_corner(p)[2])


"""
    geo_grid_centre_single(p)
    ---> Tuple{Int, Int}

When the number of sample points is even, the actual centre lies between elements.
This returns the most easterly and southerly elements in those cases.

By 'the grid' we mean internal sample points of p.

Also see `geo_centre`.
"""
geo_grid_centre_single(p) = div.(southwest_internal_corner(p) .+ northeast_external_corner(p), 2)

"""
    geo_centre(p)
    ---> Tuple{Float64, Float64}

This returns the centre point as a continuous coordinate. Since corner coordinates of
`p` are integers, the centre lies at 'n.0' for odd dimensions of `p`, or 'n.5' for even.

Also see `geo_grid_centre_single`.
"""
geo_centre(p) = (southwest_external_corner(p) .+ northeast_external_corner(p)) ./ 2

"""
    bbox_external_string(smb::SheetMatrixBuilder)
    bbox_external_string(sb::SheetBuilder)
    bbox_external_string(g::GeoArray)
    bbox_external_string(fna::String)
    ---> String

No consideration of zero padding. The large boundary box.

# Example
```
julia> bbox_external_string(smb)
"(44000 6909047)-(44018 6909063)"

julia> bbox_external_string(smb[2, 2])
"(44006 6909055)-(44012 6909063)"

julia> bbox_external_string(fna)
"(43200 6909000)-(44000 6909600)"
```
"""
function bbox_external_string(p)
    min_x, min_y = southwest_external_corner(p)
    max_x, max_y = northeast_external_corner(p)
    bbox_external_string((;min_x, min_y, max_x, max_y))
end
function bbox_external_string(bb::NamedTuple{(:min_x, :min_y, :max_x, :max_y)})
    x1 = Int(bb.min_x)
    y1 = Int(bb.min_y)
    x3 = Int(bb.max_x)
    y3 = Int(bb.max_y)
    "($x1 $y1)-($x3 $y3)"
end



"""
    geo_area(p)
    ---> Int

# Example
```
julia> geo_area(smb)
791017920

julia> northeast_external_corner(smb) .- southwest_external_corner(smb) |> tup -> *(tup...)
791017920

julia> geo_area(smb[1])
65918160
```
"""
function geo_area(p)
    s, w = southwest_external_corner(p)
    n, e = northeast_external_corner(p)
    (n - s) * (e - w)
end


"""
    closed_polygon_string(smb::SheetMatrixBuilder)
    closed_polygon_string(sb::SheetBuilder)
    closed_polygon_string(fna::String)
    closed_polygon_string(fnas::Vector{String})
    closed_polygon_string(g::GeoArray)

Used by `polygon_string`, differentiates on input type.
"""
function closed_polygon_string(smb::SheetMatrixBuilder)
    s = ""
    for sb in smb
        s *= "$(closed_polygon_string(sb))"
        if sb.sheet_number !== length(smb)
            s *= ",\n" * repeat(' ', 19)
        end
    end
    s
end
function closed_polygon_string(sb::SheetBuilder)
    min_x, min_y = southwest_external_corner(sb)
    max_x, max_y = northeast_external_corner(sb)
    closed_box_string((;min_x, min_y, max_x, max_y))
end

"""
    closed_box_string(bb::NamedTuple{(:min_x, :min_y, :max_x, :max_y)})
    ---> String

Used indirectly by `polygon_string`. Input must be convertable to integer numbers.

# Example
```
julia> bb = (;min_x = 0, min_y = 0, max_x = 1, max_y = 1);

julia> bbf = (;min_x = 0, min_y = 0, max_x = 1, max_y = 1);

julia> BitmapMaps.closed_box_string(bb)
"(0 0, 1 0, 1 1, 0 1, 0 0)"

julia> BitmapMaps.closed_box_string(bbf)
"(0 0, 1 0, 1 1, 0 1, 0 0)"
```
"""
function closed_box_string(bb::NamedTuple{(:min_x, :min_y, :max_x, :max_y)})
    x1 = Int(bb.min_x)
    y1 = Int(bb.min_y)
    x3 = Int(bb.max_x)
    y3 = Int(bb.max_y)
    x2, y2 = x3, y1
    x4, y4 = x1, y3
    "($x1 $y1, $x2 $y2, $x3 $y3, $x4 $y4, $x1 $y1)"
end


"""
    polygon_string(smb::SheetMatrixBuilder)
    polygon_string(sb::SheetBuilder)
    polygon_string(fna::String)
    polygon_string(fnas::Vector{String})
    polygon_string(g::GeoArray)
    ---> String

Paste WKS output in e.g. https://nvdb-vegdata.github.io/nvdb-visrute/STM.

# Example
```
julia> polygon_string(smb[1])
"POLYGON ((4873 6909048, 11638 6909048, 11638 6918792, 4873 6918792, 4873 6909048))"

julia> polygon_string(smb) |> println
POLYGON ((4873 6909048, 11638 6909048, 11638 6918792, 4873 6918792, 4873 6909048),
                (4873 6918792, 11638 6918792, 11638 6928536, 4873 6928536, 4873 6918792),
                (4873 6928536, 11638 6928536, 11638 6938280, 4873 6938280, 4873 6928536),
                (11638 6909048, 18403 6909048, 18403 6918792, 11638 6918792, 11638 6909048),
                (11638 6918792, 18403 6918792, 18403 6928536, 11638 6928536, 11638 6918792),
                (11638 6928536, 18403 6928536, 18403 6938280, 11638 6938280, 11638 6928536),
                (18403 6909048, 25168 6909048, 25168 6918792, 18403 6918792, 18403 6909048),
                (18403 6918792, 25168 6918792, 25168 6928536, 18403 6928536, 18403 6918792),
                (18403 6928536, 25168 6928536, 25168 6938280, 18403 6938280, 18403 6928536),
                (25168 6909048, 31933 6909048, 31933 6918792, 25168 6918792, 25168 6909048),
                (25168 6918792, 31933 6918792, 31933 6928536, 25168 6928536, 25168 6918792),
                (25168 6928536, 31933 6928536, 31933 6938280, 25168 6938280, 25168 6928536))

julia> polygon_string(fna)
"POLYGON ((44000 6909048, 44005 6909048, 44005 6909055, 44000 6909055, 44000 6909048))"
```
"""
polygon_string(p) = "POLYGON ($(closed_polygon_string(p)))"


"""
    show_augmented(smb::SheetMatrixBuilder)
    ---> Nothing

Prints to stdout, used by `run_bitmapmap_pipeline` > `define_builder`.
"""
function show_augmented(smb::SheetMatrixBuilder)
    printstyled("Bitmapmap builder based on .ini file and keywords\n", color = :green, bold=:true)
    show(stdout, MIME("text/plain"), smb)
    show_derived_properties(smb)
end

"""
    show_derived_properties(smb::SheetMatrixBuilder)
    show_derived_properties(sb::SheetBuilder)
    show_derived_properties(fna::String)
    show_derived_properties(fnas::Vector{String})
    show_derived_properties(g::GeoArray)
    ---> Nothing

The most general function for inspecting files, builders and geoarrays in this context.
Prints to stdout, used by `run_bitmapmap_pipeline` > `define_builder` > `show_augmented`.


# Example
```
julia> show_derived_properties(smb[2,2])

        Derived properties (all as (easting, northing)):
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
function show_derived_properties(p)
    tb1 = repeat(' ', 8)
    tb2 = repeat(' ', 10)
    tb3 = repeat(' ', 20)
    if p isa String
        pt = readclose(p)
    else
        pt = p
    end
    # Unit is "utm meter"
    printstyled(tb1, "\n", tb1, "Derived properties (all as (easting, northing)): \n", color = :green)
    println(tb2, rpad("Bounding Box (BB) SE-NW", 35), "= ", bbox_external_string(pt))
    println(tb2, rpad("Northeast internal corner",     35), "= ", northeast_internal_corner(pt), " - most northeastern sample point")
    println(tb2, rpad("Geo centre",        35), "= ", geo_centre(pt))
    println(tb2, rpad("Grid centre single", 35), "= ", geo_grid_centre_single(pt))
    # Other units
    printstyled(tb1, "Derived properties: \n", color = :green)
    areaint = Int(round(geo_area(pt) / 1e6))
    if areaint > 3
        println(tb2, rpad("Geographical area [km²]", 35), "= ", areaint)
    else
        println(tb2, rpad("Geographical area [m²]", 35), "= ", Int(round(geo_area(pt))))
    end
    if pt isa SheetMatrixBuilder
        a_km2 = round(geo_area(first(pt)) / 1e6, digits = 2)
        if a_km2 > 1
            println(tb3, "Per sheet = ", a_km2, " km²   (Single file export limit: 16 km²)")
        else
            println(tb3, "Per sheet = ", Int(round(geo_area(first(pt)))), " m²   (Single file export limit: 16 km²)")
        end
    end
    #
    if pt isa SheetMatrixBuilder || pt isa SheetBuilder
        w, h = round(width_adjusted_mm(pt), digits = 1), round(height_adjusted_mm(pt), digits = 1)
        println(tb2, rpad("Adj. paper (width, height) [mm]", 35), "= ", "($w, $h)")
        if pt isa SheetMatrixBuilder
            sw, sh = round(width_adjusted_mm(pt) / ncols(pt), digits = 1), round(height_adjusted_mm(pt)  / nrows(pt), digits = 1)
            println(tb3, "Per sheet [mm] (width, height) = ($sw, $sh)")
        end
    end
    #
    if pt isa SheetMatrixBuilder || pt isa SheetBuilder
        println(tb2, rpad("Map scale", 35), "= ", "1 : ", map_scale_denominator(pt), " = 1 : (cell_to_utm_factor * density_pt_m⁻¹)")
    else
        println(tb2, rpad("Distance between cells [utm or m]", 35), "= ", cell_to_utm_factor(pt))
    end
    ps = polygon_string(pt)
    if pt isa SheetMatrixBuilder
        printstyled(tb1, "BBs of sheets as Well Known Text ", color = :green)
    else
        if contains(ps, '\n')
            printstyled(tb1, "External and zero-padded boundary box as Well Known Text ", color = :green)
        else
            printstyled(tb1, "External boundary box as Well Known Text ", color = :green)
        end
    end
    printstyled("(paste in wktmap.com or nvdb-vegdata.github.io/nvdb-visrute/STM ):\n", color = :light_black)
    println(tb2, ps)
end


"""
    get_fields_namedtuple(smb::SheetMatrixBuilder)
    ---> NamedTuple

Used by `run_bitmapmap_pipeline` to deconstruct a SheetMatrixBuilder to keywords.
"""
function get_fields_namedtuple(smb::SheetMatrixBuilder)
    field_names = fieldnames(SheetMatrixBuilder)
    values = getfield.(Ref(smb), field_names)
    NamedTuple{field_names}(values)
end
