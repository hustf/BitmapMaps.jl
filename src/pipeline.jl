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
    sheet_pix_width = floor(pwi * pdens / 25.4)
    sheet_pix_height = floor(phe * pdens / 25.4)
    southwest_corner = get_config_value("Geographical position", "Southwest corner, (easting northing)", Tuple{Int, Int})
    nrows, ncols = get_config_value("Number of printable sheets", "(rows columns)", Tuple{Int, Int})
    pix_width = ncols * sheet_pix_width
    pix_height = nrows * sheet_pix_width
    pix_to_utm_factor = get_config_value("Zoom level", "Pixel distance between elevation sampling points", Int)
    bmp = BmPartition(pix_width, pix_height, sheet_pix_width, sheet_pix_height, southwest_corner, pix_to_utm_factor)
    show_augmented(bmp)
    bmp
end
