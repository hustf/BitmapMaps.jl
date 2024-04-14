"""
    run_bitmapmap(; kwds...)

The job is defined in file BitmapMaps.ini, in user's home directory.
Changing the .ini file values is recommended, although you can also 
overrule variable values using keywords. Keywords are included in
the default .ini file.

# Example
```
julia> run_bitmapmap_pipeline(nrc = (2, 2));
Bitmapmap configuration based on .ini file  BmPartition(
        pix_width                        = 4510,
        pix_height                       = 6496,
        sheet_pix_width                  = 2255,
        sheet_pix_height                 = 3248,
        nrows                            = 2,
        ncols                            = 2,
        southwest_corner                 = (4873, 6909048),
        sheet_indices                    = CartesianIndices((1:2, 1:2)),
        pix_to_utm_factor                = 3)
        Augmented properties (all as (easting, northing)):
          Geo centre =                       (11638.0, 6.918792e6)
          Grid centre single =               (11636, 6918790)
          Northeast external corner =        (18403, 6928536)
          Northeast internal corner =        (18400, 6928533) - most northeastern sample point
          Bounding Box (BB) SE-NW =          (4873 6909048)-(18403 6928536)
        BBs of sheets as Well Known Text (paste in e.g. https://nvdb-vegdata.github.io/nvdb-visrute/STM ):
          POLYGON ((4873 6909048, 11638 6909048, 11638 6918792, 4873 6918792, 4873 6909048),
                (4873 6918792, 11638 6918792, 11638 6928536, 4873 6928536, 4873 6918792),
                (11638 6909048, 18403 6909048, 18403 6918792, 11638 6918792, 11638 6909048),
                (11638 6918792, 18403 6918792, 18403 6928536, 11638 6928536, 11638 6918792))
```
"""
function run_bitmapmap_pipeline(; kwds...)
    # - Make a grid for individual printed pages. Each page is associated with an utm coordinate bounding box.
    pwi = get_kw_or_config_value(:pwi ,"Printer consistent capability", "Printable width mm", Int; kwds...)
    phe = get_kw_or_config_value(:phe ,"Printer consistent capability", "Printable height mm", Int; kwds...)
    pdensmax = get_kw_or_config_value(:pdensmax ,"Printer consistent capability", "Stated density limit, dots per inch", Int; kwds...)
    pdens = get_kw_or_config_value(:pdens ,"Printing pixel density", "Selected density, dots per inch", Int; kwds...)
    @assert pdens < pdensmax
    sheet_pix_width = floor(pwi * pdens / 25.4)
    sheet_pix_height = floor(phe * pdens / 25.4)
    southwest_corner = get_kw_or_config_value(:southwest_corner ,"Geographical position", "Southwest corner (easting northing)", Tuple{Int, Int}; kwds...)
    nrc = get_kw_or_config_value(:nrc ,"Number of printable sheets", "(rows columns)", Tuple{Int, Int}; kwds...)
    nrows, ncols = nrc
    pix_width = ncols * sheet_pix_width
    pix_height = nrows * sheet_pix_height
    pix_to_utm_factor = get_kw_or_config_value(:pix_to_utm_factor, "Pixel to utm factor", "Pixel distance between elevation sampling points", Int; kwds...)
    bmp = BmPartition(pix_width, pix_height, sheet_pix_width, sheet_pix_height, southwest_corner, pix_to_utm_factor)
    show_augmented(bmp)
    bmp
end

function 