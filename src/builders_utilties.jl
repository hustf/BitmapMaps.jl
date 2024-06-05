# Utilty functions
# Argument names 
# - smb::SheetMatrixBuilder
# - sb::SheetBuilder
# - p - duck typed, can be a GeoArray or a SheetBuilder.
#
# A few methods are extended for GeoArray in `geoarray_utilties`                 
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
northeast_corner(smb::SheetMatrixBuilder) = southwest_corner(smb) .+ smb.cell_to_utm_factor .* (smb.bm_cell_width, smb.bm_cell_height)
"ref. southwest_corner"
northwest_corner(smb::SheetMatrixBuilder) = southwest_corner(smb) .+ smb.cell_to_utm_factor .* (0, smb.bm_cell_height)
"ref. southwest_corner"
southeast_corner(smb::SheetMatrixBuilder) = southwest_corner(smb) .+ smb.cell_to_utm_factor .* (smb.bm_cell_width, 0)
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

# TODO: Restrict this method, don't accept builders.
nonzero_raster_string(p) = replace("$(southwest_external_corner(p))-$(northeast_external_corner(p))", "," => "")

function geo_area(p)
    s, w = southwest_external_corner(p)
    n, e = northeast_external_corner(p)
    (n - s) * (e - w)
end


"""
    closed_polygon_string(smb::SheetMatrixBuilder)
    closed_polygon_string(p)

Well known text (for pasting elsewhere).
"""
function closed_polygon_string(smb::SheetMatrixBuilder)
    s = ""
    for sb in smb
        s *= "$(closed_polygon_string(sb))"
        if sb.sheet_number !== length(smb)
            s *= ",\n\t\t"
        end
    end
    s
end
function closed_polygon_string(sb::SheetBuilder) 
    min_x, min_y = southwest_external_corner(sb)
    max_x, max_y = northeast_external_corner(sb)
    closed_box_string((;min_x, min_y, max_x, max_y))
end
function closed_box_string(bb::T) where T <: @NamedTuple{min_x::Float64, min_y::Float64, max_x::Float64, max_y::Float64}
    min_x = Int(bb.min_x)
    min_y = Int(bb.min_y)
    max_x = Int(bb.max_x)
    max_y = Int(bb.max_y)
    bbi = (;min_x, min_y, max_x, max_y)
    closed_box_string(bbi)
end
function closed_box_string(bb::T) where T <: @NamedTuple{min_x::Int, min_y::Int, max_x::Int, max_y::Int}
    x1 = bb.min_x
    y1 = bb.min_y
    x3 = bb.max_x
    y3 = bb.max_y
    x2, y2 = x3, y1
    x4, y4 = x1, y3
    "($x1 $y1, $x2 $y2, $x3 $y3, $x4 $y4, $x1 $y1)"
end


"""
   # polygon_string(fna::String)
    polygon_string(p)
    polygon_string(smb::SheetMatrixBuilder)

Paste output in e.g. https://nvdb-vegdata.github.io/nvdb-visrute/STM.
"""
polygon_string(p) = "POLYGON ($(closed_polygon_string(p)))"
polygon_string(fna::String) = "POLYGON ($(nonzero_raster_closed_polygon_string(fna)))"
function show_augmented(smb::SheetMatrixBuilder)
    printstyled("Bitmapmap configuration based on .ini file  ", color = :green, bold=:true)
    show(stdout, MIME("text/plain"), smb)
    show_augmented_properties(smb)
end

function show_augmented_properties(p)
    printstyled("\tAugmented properties (all as (easting, northing)): \n", color = :green)
    println("\t  ", rpad("Geo centre = ",        35), geo_centre(p))
    println("\t  ", rpad("Grid centre single = ", 35), geo_grid_centre_single(p))
    println("\t  ", rpad("Northeast external corner = ",     35), northeast_external_corner(p))
    println("\t  ", rpad("Northeast internal corner = ",     35), northeast_internal_corner(p), " - most northeastern sample point")
    println("\t  ", rpad("Bounding Box (BB) SE-NW = ", 35), nonzero_raster_string(p))
    if p isa SheetMatrixBuilder
        println("\t  ", rpad("Geographical area [km²] = ", 35), Int(round(geo_area(p) / 1e6)), "         Per sheet: ", round(geo_area(first(p)) / 1e6, digits = 1), "  km²   Single file export limit: 16 km²")
        printstyled("\tBBs of sheets as Well Known Text ", color = :green)
    else # SheetBuilder or GeoArray
        println("\t  ", rpad("Geographical area [km²] = ", 35), Int(round(geo_area(p) / 1e6)))
        printstyled("\tExternal boundary box as Well Known Text ", color = :green)
    end
    printstyled("(paste in e.g. https://nvdb-vegdata.github.io/nvdb-visrute/STM ):\n", color = :light_black)
    println("\t  ", polygon_string(p))
end

function get_fields_namedtuple(smb::SheetMatrixBuilder)
    field_names = fieldnames(SheetMatrixBuilder)
    values = getfield.(Ref(smb), field_names)
    NamedTuple{field_names}(values)
end