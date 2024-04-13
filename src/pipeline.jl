"""
    run_bitmapmap()

The job is defined in file BitmapMaps.ini, in user's home directory.

"""
function run_bitmapmap_pipeline()
    # - Make a grid for individual printed pages. Each page is associated with an utm coordinate bounding box.
    pwi = get_config_value("Printer consistent capability", "Printable width mm", Int)
    phe = get_config_value("Printer consistent capability", "Printable height mm", Int)
    pdensmax = get_config_value("Printer consistent capability", "Stated density limit, dots per inch", Int)
    pdens = get_config_value("Printing pixel density", "Selected density, dots per inch", Int)
    @assert pdens < pdensmax
    sheet_width = floor(pwi * pdens / 25.4)
    sheet_height = floor(phe * pdens / 25.4)
    bm_southwest_corner = get_config_value("Geographical position", "Southwest corner, (easting northing)", Tuple{Int, Int})
    nrows, ncols = get_config_value("Number of printable sheets", "(rows columns)", Tuple{Int, Int})
    bm_pixel_width = ncols * sheet_width
    bm_pixel_height = nrows * sheet_width
    pixel_distance = get_config_value("Zoom level", "Pixel distance between elevation sampling points", Int)
    bmp = BmPartition(bm_pixel_width, bm_pixel_height, sheet_width, sheet_height, bm_southwest_corner, pixel_distance)
    show_augmented(bmp)
    bmp
end
function show_augmented(p::BmPartition)
    printstyled("Bitmapmap configuration based on .ini file  ", color = :green, bold=:true)
    show(stdout, MIME("text/plain"), p)
    printstyled("\tAugmented properties (all as (easting, northing)): \n", color = :green)
    println("\t  ", rpad("Geo centre = ",        30), geo_centre(p))
    println("\t  ", rpad("Grid centre single = ", 30), geo_grid_centre_single(p))
    println("\t  ", rpad("Northeast external corner = ",        30), northeast_external_corner(p))
    println("\t  ", rpad("Northeast internal corner = ",        30), northeast_internal_corner(p), " - most northeastern sample point")
    println("\t  ", rpad("Bounding box SE-NW (external) = ",        30), bounding_box_external_string(p))
end