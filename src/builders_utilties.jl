# Utilty functions
"""
    southwest_corner(p::SheetMatrixBuilder)
    northeast_corner(p::SheetMatrixBuilder)
    northwest_corner(p::SheetMatrixBuilder) 
    southeast_corner(p::SheetMatrixBuilder) 
    northwest_corner(p::SheetBuilder)
    southeast_internal_corner(p::SheetBuilder)
    southeast_external_corner(p::SheetBuilder)
    southwest_internal_corner(p::SheetBuilder)
    southwest_external_corner(p::SheetBuilder)
    northeast_internal_corner(p::SheetBuilder)
    northeast_external_corner(p::SheetBuilder)
    ---> Tuple(Int, Int)

An entire SheetMatrixBuilder is divided into sheets. Every SheetMatrixBuilder method refers to 'external' borders.
A SheetBuilder is a subset of a BmParition. It shares corners with neighbouring sheets.
Each sheet is divided into pixels, and each pixel refers to a single geographical position or sampling point.
    
'external' refers to outside dimensions of the sheet or BitMap, which are larger than 'internal'.
'internal' refers to the last sampling point within or on the border of 'external' borders.

Each geographical sampling point corresponds to the top left of its pixel. Hence, 'external' and 'internal'
are indentical for the northwest corner of a sheet.
"""
southwest_corner(p::SheetMatrixBuilder) = p.southwest_corner
"ref. southwest_corner"
northeast_corner(p::SheetMatrixBuilder) = southwest_corner(p) .+ p.pix_to_utm_factor .* (p.bm_pix_width, p.bm_pix_height)
"ref. southwest_corner"
northwest_corner(p::SheetMatrixBuilder) = southwest_corner(p) .+ p.pix_to_utm_factor .* (0, p.bm_pix_height)
"ref. southwest_corner"
southeast_corner(p::SheetMatrixBuilder) = southwest_corner(p) .+ p.pix_to_utm_factor .* (p.bm_pix_width, 0)
"ref. southwest_corner"
northwest_corner(p::SheetBuilder) = p.f_I_to_utm(CartesianIndex(1, 1))
southwest_external_corner(p::SheetMatrixBuilder) = southwest_corner(p)
northeast_external_corner(p::SheetMatrixBuilder) = northeast_corner(p)
northwest_external_corner(p::SheetMatrixBuilder) = northwest_corner(p)
southeast_external_corner(p::SheetMatrixBuilder) = southeast_corner(p)
"ref. southwest_corner"
southeast_internal_corner(p::SheetBuilder) = p.f_I_to_utm(p.pix_iter[end, end])
southeast_internal_corner(p::SheetMatrixBuilder) = southeast_internal_corner(p[1, end])
"ref. southwest_corner"
southeast_external_corner(p::SheetBuilder) = 2 .* southeast_internal_corner(p) .- p.f_I_to_utm(p.pix_iter[end - 1, end - 1])
"ref. southwest_corner"
southwest_internal_corner(p::SheetBuilder) = (northwest_corner(p)[1], southeast_internal_corner(p)[2])
southwest_internal_corner(p::SheetMatrixBuilder) = southwest_internal_corner(p[1, 1])
"ref. southwest_corner"
southwest_external_corner(p::SheetBuilder) = (northwest_corner(p)[1], southeast_external_corner(p)[2])
"ref. southwest_corner"
northeast_internal_corner(p::SheetBuilder) = (southeast_internal_corner(p)[1], northwest_corner(p)[2])
northeast_internal_corner(p::SheetMatrixBuilder) = northeast_internal_corner(p[end, end])
"ref. southwest_corner"
northeast_external_corner(p::SheetBuilder) = (southeast_external_corner(p)[1], northwest_corner(p)[2])


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

bounding_box_external_string(p) = replace("$(southwest_external_corner(p))-$(northeast_external_corner(p))", "," => "")

function geo_area(p)
    s, w = southwest_external_corner(p)
    n, e = northeast_external_corner(p)
    (n - s) * (e - w)
end

# Well known text (for pasting elsewhere)

function bounding_box_closed_polygon_string(p::T) where T <: Union{SheetBuilder, SheetMatrixBuilder}
    x1, y1 = southwest_external_corner(p)
    x3, y3 = northeast_external_corner(p)
    x2, y2 = x3, y1
    x4, y4 = x1, y3
   "($x1 $y1, $x2 $y2, $x3 $y3, $x4 $y4, $x1 $y1)"
end


function bounding_box_closed_polygon_string(smb::SheetMatrixBuilder)
    s = ""
    for sb in smb
        s *= "$(bounding_box_closed_polygon_string(sb))"
        if sb.sheet_number !== length(smb)
            s *= ",\n\t\t"
        end
    end
    s
end

function bounding_box_closed_polygon_string(fna::String)
    @assert isfile(fna)
    @assert endswith(fna, r".tif|.TIF")
    bounding_box_closed_polygon_string(bbox(GeoArrays.read(fna)))
end
function bounding_box_closed_polygon_string(bb::T) where T <: @NamedTuple{min_x::Float64, min_y::Float64, max_x::Float64, max_y::Float64}
    x1 = Int(bb.min_x)
    y1 = Int(bb.min_y)
    x3 = Int(bb.max_x)
    y3 = Int(bb.max_y)
    x2, y2 = x3, y1
    x4, y4 = x1, y3
   "($x1 $y1, $x2 $y2, $x3 $y3, $x4 $y4, $x1 $y1)"
end




"""
    bounding_box_polygon_string(fna::String)
    bounding_box_polygon_string(p::SheetBuilder)
    bounding_box_polygon_string(p::SheetMatrixBuilder)

Paste output in e.g. https://nvdb-vegdata.github.io/nvdb-visrute/STM.
"""
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

