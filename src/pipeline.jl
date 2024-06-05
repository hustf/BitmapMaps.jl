"""
    run_bitmapmap_pipeline(; complete_sheets_first = true, kwds...)
    ---> SheetMatrixBuilder

The job is defined in file BitmapMaps.ini, in user's home directory.
You can overrule parameters with keywords, but changing the .ini file 
is recommended.

# Arguments

`complete_sheets_first`    The default 'true' means one sheet is fully processed, then the next sheet.
                         'false' means that each operation is finished for all sheets before the next operation. 
`bm_cell_width`, etc           Keyword names (like `bm_cell_width`) are included in the default .ini file, with explanation.

# Example
```
julia> run_bitmapmap_pipeline(;nrc = (2, 2));
Bitmapmap configuration based on .ini file  SheetMatrixBuilder(     4510, # bm_cell_width
                        6496, # bm_cell_height
                        2255, # sheet_cell_width
                        3248, # sheet_cell_height
                           2, # nrows
                           2, # ncols
             (4873, 6909048), # southwest_corner
CartesianIndices((1:2, 1:2)), # sheet_indices
                           3, # cell_to_utm_factor
        "bitmapmaps/default") # pth
        Augmented properties (all as (easting, northing)):
          Geo centre =                       (11638.0, 6.918792e6)
          Grid centre single =               (11638, 6918793)
          Northeast external corner =        (18403, 6928536)
          Northeast internal corner =        (18400, 6928536) - most northeastern sample point
          Bounding Box (BB) SE-NW =          (4873 6909048)-(18403 6928536)
          Geographical area [km²] =          264         Per sheet: 65.9  km²   Single file export limit: 16 km²
        BBs of sheets as Well Known Text (paste in e.g. https://nvdb-vegdata.github.io/nvdb-visrute/STM ):
          POLYGON ((4873 6909048, 11638 6909048, 11638 6918792, 4873 6918792, 4873 6909048),
                (4873 6918792, 11638 6918792, 11638 6928536, 4873 6928536, 4873 6918792),
                (11638 6909048, 18403 6909048, 18403 6918792, 11638 6918792, 11638 6909048),
                (11638 6918792, 18403 6918792, 18403 6928536, 11638 6928536, 11638 6918792))
```
"""
function run_bitmapmap_pipeline(; complete_sheets_first = true, kwds...)
    smb = define_builder(; kwds...)
    ok_res = process_job(smb, complete_sheets_first)
    if ok_res
        printstyled("Finished job in folder $(joinpath(homedir(), smb.pth))\n", color = :yellow)
    end
    smb
end
function run_bitmapmap_pipeline(smb::SheetMatrixBuilder; kwds...)
    # Extract current fields as a NamedTuple.
    current_fields = get_fields_namedtuple(smb)
    # Combine current fields with new keyword arguments.
    # New keyword arguments would overrule.
    combined_kwds = merge(current_fields, kwds)
    # Call the original define_builder with the combined keyword arguments.
    run_bitmapmap_pipeline(; combined_kwds...)
end



function define_builder(; kwds...)
    # - Make a grid for individual printed pages. Each page is associated with an utm coordinate bounding box.
    pwi = get_kw_or_config_value(:pwi ,"Printer consistent capability", "Printable width mm", Int; kwds...)
    phe = get_kw_or_config_value(:phe ,"Printer consistent capability", "Printable height mm", Int; kwds...)
    pdensmax_dpi = get_kw_or_config_value(:pdensmax_dpi ,"Printer consistent capability", "Stated density limit, dots per inch", Int; kwds...)
    pdens_dpi = get_kw_or_config_value(:pdens_dpi ,"Printing pixel density", "Selected density, dots per inch", Int; kwds...)
    southwest_corner = get_kw_or_config_value(:southwest_corner ,"Geographical position", "Southwest corner (easting northing)", Tuple{Int, Int}; kwds...)
    nrc = get_kw_or_config_value(:nrc ,"Number of printable sheets", "(rows columns)", Tuple{Int, Int}; kwds...)
    cell_to_utm_factor = get_kw_or_config_value(:cell_to_utm_factor, "Cell to utm factor", "Utm unit distance between elevation sampling points", Int; kwds...)
    @assert pdens_dpi < pdensmax_dpi
    pth = get_kw_or_config_value(:pth, "File folder", "Top folders path under homedir()", String; kwds...)
    if isnothing(match(r"(b|B)itmap(M|m)aps", pth))
        throw(ArgumentError("Any BmParition path must match regex r\"(b|B)itmap(M|m)aps\". Current path is $pth"))
    end
    smb = SheetMatrixBuilder(pwi, phe, pdens_dpi, nrc, southwest_corner, cell_to_utm_factor, pth)
    show_augmented(smb)
    if geo_area(smb) >= 16e6 
        @info "Since geographical area per sheet is  > 16km², 'høydedata.no' will not provide a single elevation data file per sheet. The pipieline will make a consolidated single file."
    end
    smb
end


function process_job(smb, complete_sheets_first)
    #= TODO: Good func names for:
    Make topographic reliefs
    Make vector graphics and text covering the full map area. You may use RouteMap.jl for this step.
    Make composite bitmaps
    =#
    operations_order = [establish_folder, unzip_tif, consolidate_elevation_data, water_overlay, topo_relief, contour_lines_overlay]
    # Consider sharing an in-memory z-map and gradient map between operations.... No need to redo it.
    if complete_sheets_first
        for sb in smb # Might do this in parallel? Though a lot will be wating for file i/o...
            for fn in operations_order
                call_func(fn, sb) || return false
            end
        end
    else
        for fn in operations_order
            for sb in smb
                call_func(fn, sb) || return false
            end
        end
    end
    true
end
function call_func(fn, sb)
    ok_res = fn(sb)
    if ! ok_res
        @warn "Could not finish $fn($sb) with success. Exiting."
        return false
    else
        @debug "Finished $fn"
    end
    true
end