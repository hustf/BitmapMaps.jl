# Utilty functions
# Argument names
# - smb::SheetMatrixBuilder
# - sb::SheetBuilder
# - p - duck typed, can be a GeoArray or a SheetBuilder.
#
# A few functions are extended with methods for GeoArray in `geoarray_utilties`


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

cell_to_utm_factor(smb::SheetMatrixBuilder) = smb.cell_to_utm_factor
cell_to_utm_factor(sb::SheetBuilder) = cell_to_utm_factor(sb.f_I_to_utm)
cell_to_utm_factor(f_I_to_utm::Function) = f_I_to_utm(CartesianIndex(1, 2))[1] - f_I_to_utm(CartesianIndex(1, 1))[1]
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
northwest_corner(sb::SheetBuilder)
southwest_external_corner(smb::SheetMatrixBuilder)
northeast_external_corner(smb::SheetMatrixBuilder)
northwest_external_corner(smb::SheetMatrixBuilder)
southeast_external_corner(smb::SheetMatrixBuilder)
southeast_internal_corner(sb::SheetBuilder)
southeast_internal_corner(smb::SheetMatrixBuilder)
southeast_external_corner(sb::SheetBuilder)
southwest_internal_corner(p)
southwest_internal_corner(smb::SheetMatrixBuilder)
southwest_external_corner(p)
northeast_internal_corner(p)
northeast_internal_corner(smb::SheetMatrixBuilder)
northeast_external_corner(p)
    ---> Tuple(Int, Int)

An entire SheetMatrixBuilder is divided into sheets. Every SheetMatrixBuilder method refers to 'external' boundaries.
A SheetBuilder is a subset of a SheetMatrixBuilder. It shares corners with neighbouring sheets.
Each sheet is divided into cells (i.e. pixels), and each cell is represented by a single geographical position or sampling point.

'external' refers to outside dimensions of the sheet or BitMap, which are larger than 'internal'.
'internal' refers to the last sampling point within or on the 'external' boundaries.

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
southeast_internal_corner(sb::SheetBuilder) = sb.f_I_to_utm(sb.cell_iter[end, end])
"ref. southwest_corner"
southeast_internal_corner(smb::SheetMatrixBuilder) = southeast_internal_corner(smb[1, end])
"ref. southwest_corner"
southeast_external_corner(sb::SheetBuilder) = 2 .* southeast_internal_corner(sb) .- sb.f_I_to_utm(sb.cell_iter[end - 1, end - 1])
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
bbox_external_string(p) = _bbox_external_string(bbox_internal(p))
function _bbox_external_string(bb::NamedTuple{(:min_x, :min_y, :max_x, :max_y)})
    x1 = Int(bb.min_x)
    y1 = Int(bb.min_y)
    x3 = Int(bb.max_x)
    y3 = Int(bb.max_y)
    "($x1 $y1)-($x3 $y3)"
end

"""
    bbox_external(p)
    --> NamedTuple{(:min_x, :min_y, :max_x, :max_y)}

UTM coordinates from p, a SheetMatrixBuilder or SheetBuilder.
"""
function bbox_external(p)
    min_x, min_y = southwest_external_corner(p)
    max_x, max_y = northeast_external_corner(p)
    (;min_x, min_y, max_x, max_y)
end
"""
    bbox_internal(p)
    bbox_internal(cell_iter, f_I_to_utm)
    --> NamedTuple{(:min_x, :min_y, :max_x, :max_y)}

UTM coordinates from p, a SheetMatrixBuilder or SheetBuilder.
"""
function bbox_internal(p)
    min_x, min_y = southwest_internal_corner(p)
    max_x, max_y = northeast_internal_corner(p)
    (;min_x, min_y, max_x, max_y)
end
function bbox_internal(cell_iter, f_I_to_utm)
    nw = f_I_to_utm(CartesianIndex(1, 1))
    se = f_I_to_utm(CartesianIndices(cell_iter)[end])
    min_x, min_y, max_x, max_y = nw[1], se[2], se[1], nw[2]
    (;min_x, min_y, max_x, max_y)
end
function bbox_internal(v::Vector)
    vbb = bbox_internal.(v)
    max_x = maximum(nt -> nt.max_x, vbb)
    min_x = minimum(nt -> nt.min_x, vbb)
    max_y = maximum(nt -> nt.max_y, vbb)
    min_y = minimum(nt -> nt.min_y, vbb)
    (;min_x, min_y, max_x, max_y)
end
function bbox_internal(g::MetaGraph)
    isempty(g[]) && throw(ArgumentError("No metadata"))
    g[] isa Set{@NamedTuple{min_x::Int64, min_y::Int64, max_x::Int64, max_y::Int64}} || throw(ArgumentError("Metadata not recognized"))
    max_x = maximum(nt -> nt.max_x, g[])
    min_x = minimum(nt -> nt.min_x, g[])
    max_y = maximum(nt -> nt.max_y, g[])
    min_y = minimum(nt -> nt.min_y, g[])
    (;min_x, min_y, max_x, max_y)
end
function bbox_internal(vtup::Vector{Tuple{Int64, Int64}})
    max_x = maximum(tup -> tup[1], vtup)
    min_x = minimum(tup -> tup[1], vtup)
    max_y = maximum(tup -> tup[2], vtup)
    min_y = minimum(tup -> tup[2], vtup)
    (;min_x, min_y, max_x, max_y)
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
    geo_width_height(p)
    ---> Int

# Example
```
julia> geo_width_height(smb)
(27060, 29232)

julia> geo_width_height(smb[1,1])
(6765, 9744)
```
"""
function geo_width_height(p)
    s, w = southwest_external_corner(p)
    n, e = northeast_external_corner(p)
    (n - s, e - w)
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
"((0 0, 1 0, 1 1, 0 1, 0 0))"

julia> BitmapMaps.closed_box_string(bbf)
"((0 0, 1 0, 1 1, 0 1, 0 0))"
```
"""
function closed_box_string(bb::NamedTuple{(:min_x, :min_y, :max_x, :max_y)})
    x1 = Int(bb.min_x)
    y1 = Int(bb.min_y)
    x3 = Int(bb.max_x)
    y3 = Int(bb.max_y)
    x2, y2 = x3, y1
    x4, y4 = x1, y3
    "(($x1 $y1, $x2 $y2, $x3 $y3, $x4 $y4, $x1 $y1))"
end

"""
    diagonal_string(bb::NamedTuple{(:min_x, :min_y, :max_x, :max_y)})
    ---> String

# Example
```
julia> diagonal_string((; min_x = 0, min_y = 0, max_x = 100, max_y = 200))
"((0 0, 1 0, 99 200, 100 200, 0 0))"
```
"""
function diagonal_string(bb::NamedTuple{(:min_x, :min_y, :max_x, :max_y)})
    x1 = Int(bb.min_x)
    y1 = Int(bb.min_y)
    x3 = Int(bb.max_x) - 1
    y3 = Int(bb.max_y)
    x2 = x1 + 1
    y2 = y1
    x4 = x3 + 1
    y4 = y3
    "(($x1 $y1, $x2 $y2, $x3 $y3, $x4 $y4, $x1 $y1))"
end


"""
    polygon_string(smb::SheetMatrixBuilder)
    polygon_string(sb::SheetBuilder)
    polygon_string(fna::String)
    polygon_string(fnas::Vector{String})
    polygon_string(g::GeoArray)
    ---> String

Paste WKS output in e.g. https://nvdb- ->ZXCVBNM;:_.github.io/nvdb-visrute/STM.

# Example
```
julia> polygon_string(smb[1])
"MULTIPOLYGON (\n                   ((35425 6920995, 42190 6920995, 42190 6930739, 35425 6930739, 35425 6920995)))"

julia> polygon_string(smb) |> println
MULTIPOLYGON (
                   ((35425 6920995, 42190 6920995, 42190 6930739, 35425 6930739, 35425 6920995)),
                   ((35425 6930739, 42190 6930739, 42190 6940483, 35425 6940483, 35425 6930739)),
                   ((35425 6940483, 42190 6940483, 42190 6950227, 35425 6950227, 35425 6940483)),
                   ((42190 6920995, 48955 6920995, 48955 6930739, 42190 6930739, 42190 6920995)),
                   ((42190 6930739, 48955 6930739, 48955 6940483, 42190 6940483, 42190 6930739)),
                   ((42190 6940483, 48955 6940483, 48955 6950227, 42190 6950227, 42190 6940483)),
                   ((48955 6920995, 55720 6920995, 55720 6930739, 48955 6930739, 48955 6920995)),
                   ((48955 6930739, 55720 6930739, 55720 6940483, 48955 6940483, 48955 6930739)),
                   ((48955 6940483, 55720 6940483, 55720 6950227, 48955 6950227, 48955 6940483)),
                   ((55720 6920995, 62485 6920995, 62485 6930739, 55720 6930739, 55720 6920995)),
                   ((55720 6930739, 62485 6930739, 62485 6940483, 55720 6940483, 55720 6930739)),
                   ((55720 6940483, 62485 6940483, 62485 6950227, 55720 6950227, 55720 6940483)))

julia> polygon_string(fna) |> println  # This file contains only zero values, which is indicated by a diagonal:
MULTIPOLYGON (
                   ((-54575 6890995, -39565 6890995, -39565 6906005, -54575 6906005, -54575 6890995)),
                   ((-54575 6890995, -54574 6890995, -39566 6906005, -39565 6906005, -54575 6890995)))
```
"""
polygon_string(p) = "MULTIPOLYGON (\n" * repeat(' ', 18) * " $(closed_polygon_string(p)))"


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
        
        [easting, northing] derived properties:
          Bounding Box (BB) SE-NW            = (42190 6930739)-(48955 6940483)
          Northeast internal corner          = (48952, 6940483) - most northeastern sample point
          Geo centre                         = (45572.5, 6.935611e6)
          Grid centre single                 = (45572, 6935612)
        Derived properties:
          Geographical (width, height) [km]  = (6.8, 9.7)
          Geographical area [km²]            = 66
          Sheets total (width, height) [cm]  = (19.1, 27.5)
          Map scale                          = 1 : 35433 = 1 : (cell_to_utm_factor * density_pt_m⁻¹)
        External and zero-padded boundary box as Well Known Text (paste in wktmap.com or nvdb-vegdata.github.io/nvdb-visrute/STM ):
          MULTIPOLYGON (
                   ((42190 6930739, 48955 6930739, 48955 6940483, 42190 6940483, 42190 6930739)))

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
    printstyled(tb1, "\n", tb1, "[easting, northing] derived properties: \n", color = :green)
    println(tb2, rpad("Bounding Box (BB) SE-NW", 35), "= ", bbox_external_string(pt))
    println(tb2, rpad("Northeast internal corner",     35), "= ", northeast_internal_corner(pt), " - most northeastern sample point")
    println(tb2, rpad("Geo centre",        35), "= ", geo_centre(pt))
    println(tb2, rpad("Grid centre single", 35), "= ", geo_grid_centre_single(pt))
    # Other units
    printstyled(tb1, "Derived properties: \n", color = :green)
    areaint = Int(round(geo_area(pt) / 1e6))
    if areaint > 3
        println(tb2, rpad("Geographical (width, height) [km]", 35), "= ", round.(geo_width_height(pt) ./ 1000, digits = 1))
        println(tb2, rpad("Geographical area [km²]", 35), "= ", areaint)
    else
        println(tb2, rpad("Geographical (width, height) [m]", 35), "= ", geo_width_height(pt))
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
        w_mm, h_mm = round(width_adjusted_mm(pt), digits = 1), round(height_adjusted_mm(pt), digits = 1)
        w_cm, h_cm = round(w_mm / 10, digits = 1), round(h_mm / 10, digits = 1)
        println(tb2, rpad("Sheets total (width, height) [cm]", 35), "= ", "($(w_cm), $(h_cm))")
        if pt isa SheetMatrixBuilder
            sw, sh = round(width_adjusted_mm(pt) / ncols(pt), digits = 1), round(height_adjusted_mm(pt)  / nrows(pt), digits = 1)
            println(tb3, "Per sheet [mm] (w, h) = ($sw, $sh)")
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


"""
    cartesian_index_string(smb::SheetMatrixBuilder, i::Int)
        cartesian_index_string(smb::SheetMatrixBuilder)
    ---> String

Used for feedback by pipeline.

# Example
```
julia> BitmapMaps.cartesian_index_string(smb)
"[7, 7]"

julia> BitmapMaps.cartesian_index_string(smb,2)
"[2, 1]"
```
"""
cartesian_index_string(smb::SheetMatrixBuilder, i::Int) = replace(string(CartesianIndices(axes(smb))[i].I), '(' => '[', ')' => ']')
cartesian_index_string(smb::SheetMatrixBuilder) = replace(string(CartesianIndices(axes(smb))[end].I), '(' => '[', ')' => ']')


"""
    func_utm_to_sheet_index(smb::SheetMatrixBuilder)

Generate a function which returns on which sheet the utm coordinate belongs. Also see
`func_I_to_utm` and `func_utm_to_cell_index`.
"""
function func_utm_to_sheet_index(smb::SheetMatrixBuilder)
    sw_easting, sw_northing = southwest_external_corner(smb)
    nr = nrows(smb)
    nc = ncols(smb)
    ws_utm, hs_utm = geo_width_height(smb[1]) # Sheet width, height
    f = let sw_easting = sw_easting, sw_northing = sw_northing, nr = nr,  nc = nc, hs_utm = hs_utm, ws_utm = ws_utm
        utm::Tuple{Int64, Int64} -> let 
            i, j = (div((utm[2] - sw_northing), hs_utm, RoundUp), div((utm[1] - sw_easting), ws_utm) + 1)
            i < 1 && throw(ArgumentError("Northing $(utm[2]) would coorespond to sheet [$i, j]. sw_northing = $(sw_northing) and hs_utm = $(hs_utm)"))
            j < 1 && throw(ArgumentError("Easting $(utm[1]) would coorespond to sheet [i, $j]. sw_easting = $(sw_easting) and ws_utm = $(hs_utm)"))
            i > nr && throw(ArgumentError("Northing $(utm[2]) would coorespond to sheet [$i, j], while nr = $nr. sw_northing = $(sw_northing) and hs_utm = $(hs_utm)"))
            j > nc && throw(ArgumentError("Easting $(utm[1]) would coorespond to sheet [i, $j], while nc = $nc. sw_easting = $(sw_easting) and ws_utm = $(hs_utm)"))
            return (i, j)
        end
    end
    f
end

"""
    func_utm_to_cell_index(p)

p can be a SheetMatrixBuilder or a SheetBuilder.

Generate a function which returns the sheet local index corresponding to the utm coordinate. 
Also see `func_I_to_utm` and `func_utm_to_sheet_index`.
"""
function func_utm_to_cell_index(p)
    # In the northwest, external and internal corners are the same.
    nw_easting, nw_northing = northwest_corner(p)
    if p isa SheetMatrixBuilder
        w_utm, h_utm = geo_width_height(p[1]) # Sheet width, height
    else
        w_utm, h_utm = geo_width_height(p) # Sheet width, height
    end
    cell2utm = cell_to_utm_factor(p)
    f = let nw_easting = nw_easting, nw_northing = nw_northing, w_utm = w_utm, h_utm = h_utm, cell2utm = cell2utm
        utm::Tuple{Int64, Int64} -> let 
            i_utm, j_utm = (mod(nw_northing - utm[2], h_utm) + 1, mod(utm[1] - nw_easting, w_utm) + 1)
            return div(i_utm, cell2utm, RoundUp), div(j_utm, cell2utm, RoundUp)
        end
    end
    f
end



# NOTE we have kept 'neighbor_folder_dict'.
# We may perhaps reuse neighbor_folder if we will be linking sheet svgs to each other. 
# For other uses, see `sides_with_border(smb, sb)`

"""
    neighbor_folder(sb::SheetBuilder, direction::Symbol)
    ---> String

Returns an empty string if the neigbour folder does not exist.
"""
function neighbor_folder(sb::SheetBuilder, direction::Symbol)
    # Prepare
    Δi, Δj = Tuple(cartesian_unit_offset(direction))
    # Collect from sb
    fo = joinpath(splitpath(sb.pthsh)[1:end - 1])
    i = parse(Int, split(splitpath(sb.pthsh)[end])[1])
    j = parse(Int, split(splitpath(sb.pthsh)[end])[2])
    sw = southwest_external_corner(sb)
    ne = northeast_external_corner(sb)
    w = width_cell(sb) * cell_to_utm_factor(sb)
    h = height_cell(sb) * cell_to_utm_factor(sb)
    #
    fofo = ""
    fofo *= string(i + Δi) * " " * string(j + Δj)
    fofo *= "  " * string(sw[1] + Δj * w) * " " * string(sw[2] + Δi * h)
    fofo *= "  " * string(ne[1] + Δj * w) * " " * string(ne[2] + Δi * h)
    fullfo = joinpath(homedir(), fo, fofo)
    if isdir(fullfo)
        return fullfo
    else
        return ""
    end
end
function neighbor_folder_dict(sb)
    dic = Dict{Symbol, String}()
    for direction in [:s, :n, :w, :e]
        nefo = neighbor_folder(sb, direction)
        if nefo !== ""
            push!(dic, direction => nefo)
        end
    end
    dic
end


"""
    sides_with_border(smb, sb)
    sides_with_border(smb, i::Int64)
    --> Vector{Symbol}

# Example
```
julia> sides_with_border(smb, 1)
2-element Vector{Symbol}:
 :n
 :e

julia> sides_with_border(smb, 9)
2-element Vector{Symbol}:
 :w
 :s

julia> sides_with_border(smb, smb[2,2])
4-element Vector{Symbol}:
 :w
 :n
 :s
 :e
```
"""
function sides_with_border(smb, i::Int64)
    nr = nrows(smb)
    nc = ncols(smb)
    borders = Symbol[]
    R = CartesianIndices((nr, nc))
    I = R[i]
    for direction in [:w, :n, :s, :e]
        Δ = cartesian_unit_offset(direction)
        I_neighbor = I + Δ
        if I_neighbor ∈ R
            push!(borders, direction)
        end
    end
    borders
end
sides_with_border(smb, sb) = sides_with_border(smb, sb.sheet_number)

function cartesian_unit_offset(direction)
    if direction == :n
        CartesianIndex(1, 0)
    elseif direction == :s
        CartesianIndex(-1, 0)
    elseif direction == :e
        CartesianIndex(0, 1)
    elseif direction == :w
        CartesianIndex(0, -1)
    else
        throw(ArgumentError("Direction is $direction, ∉ [:n, :s, :e, :w]"))
    end
end