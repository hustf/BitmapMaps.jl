# This pipeline is the most central part of the package.
"""
    run_bitmapmap_pipeline(; complete_sheets_first = true, kwds...)
    ---> SheetMatrixBuilder

The job is defined in file BitmapMaps.ini, in user's home directory.
You can overrule parameters from the .ini file with identical keywords, but changing the .ini file
is recommended.

# Arguments

`complete_sheets_first`    The default 'true' means one sheet is fully processed, then the next sheet.
                         'false' means that each operation is finished for all sheets before the next operation.
`width_cell`, etc           Keyword names (like `width_cell`) are included in the default .ini file, with explanation.

# Example
```
julia> run_bitmapmap_pipeline(;nrc = (2, 2));
Bitmapmap builder based on .ini file and keywords
SheetMatrixBuilder((4873, 6909048), # southwest_corner
      CartesianIndices((1:2, 1:2)), # sheet_indices
                                 3, # cell_to_utm_factor
                               191, # sheet_width_mm
                               275, # sheet_height_mm
                             11811, # density_pt_m⁻¹
              "bitmapmaps/default") # pth
        [easting, northing] derived properties:
          Bounding Box (BB) SE-NW            = (4873 6909048)-(18403 6928536)
          Northeast internal corner          = (18400, 6928536) - most northeastern sample point
          Geo centre                         = (11638.0, 6.918792e6)
          Grid centre single                 = (11638, 6918793)
        Derived properties:
          Geographical area [km²]            = 264
                    Per sheet = 65.92 km²   (Single file export limit: 16 km²)
          Adj. paper (width, height) [mm]    = (381.8, 550.0)
                    Per sheet [mm] (width, height) = (190.9, 275.0)
          Map scale                          = 1 : 35433 = 1 : (cell_to_utm_factor * density_pt_m⁻¹)
        BBs of sheets as Well Known Text (paste in e.g. https://nvdb-vegdata.github.io/nvdb-visrute/STM ):
          POLYGON ((4873 6909048, 11638 6909048, 11638 6918792, 4873 6918792, 4873 6909048),
                   (4873 6918792, 11638 6918792, 11638 6928536, 4873 6928536, 4873 6918792),
                   (11638 6909048, 18403 6909048, 18403 6918792, 11638 6918792, 11638 6909048),
                   (11638 6918792, 18403 6918792, 18403 6928536, 11638 6928536, 11638 6918792))
[ Info: Since geographical area per sheet is  > 16km², 'høydedata.no' will not provide a single elevation data file per sheet. The pipeline will make a consolidated single file.
[ Info: No relevant data to consolidate. Exiting.
[ Info: Could not make consolidated .tif for sheet  with folder path bitmapmaps/default\1 1  4873 6909048  11638 6918792. Download and unzip .tif files? Exiting.
┌ Warning: Could not finish consolidate_elevation_data(SheetBuilder((0, 3248), (1:3248, 1:2255), f(I) -> utm, 111811, "bitmapmaps/default\\1 1  4873 6909048  11638 6918792",  )
│ ) with success. Exiting.
```
"""
function run_bitmapmap_pipeline(; complete_sheets_first = true, kwds...)
    if get(ENV, "JULIA_DEBUG", "") !== "BitmapMaps"
        @info "Pipeline running. Set ENV[\"JULIA_DEBUG\"] = \"BitmapMaps\" for detailed progress."
    end
    smb = define_builder(; kwds...)
    ok_res = process_job(smb, complete_sheets_first)
    if ok_res
        printstyled("Finished job in folder $(joinpath(homedir(), smb.pth))\n", color = :yellow)
        for sb in smb
            printstyled(repeat(' ', 56), sb.pthsh[length(smb.pth) + 1:end], "\n", color = :yellow)
        end
    end
    smb
end
function run_bitmapmap_pipeline(smb::SheetMatrixBuilder; kwds...)
    # Deconstruct SheetMatrixBuilder to keywords.
    current_fields = get_fields_namedtuple(smb)
    # Combine current fields with new keyword arguments.
    # New keyword arguments would overrule.
    combined_kwds = merge(current_fields, kwds)
    # Run pipeline with the combined keyword arguments.
    run_bitmapmap_pipeline(; combined_kwds...)
end


"""
    function define_builder(; kwds...)

Make a grid for individual printed pages. Each page is associated with an utm coordinate bounding box.

Parameters from .ini file will be overridden by keywords. See `run_bitmapmap_pipeline`.
"""
function define_builder(; kwds...)
    allowed_keywords = [:southwest_corner, :cell_to_utm_factor, :sheet_width_mm, :sheet_height_mm, :density_pt_m⁻¹, :pth, :density_limit_pt_inch⁻¹, :sheet_indices, :nrc]
    unrecognized_keywords = filter(∉(allowed_keywords), keys(kwds))
    if ! isempty(unrecognized_keywords)
        throw(ArgumentError("Unrecognized_keywords: $unrecognized_keywords. See file BitmapMap.ini for keywords, like: 'sheet_width_mm'"))
    end
    # Parameters from .ini file, overridden by key words.
    southwest_corner = get_kw_or_config_value(:southwest_corner ,"Geographical area", "Southwest corner (utm easting northing)", Tuple{Int, Int}; kwds...)
    cell_to_utm_factor = get_kw_or_config_value(:cell_to_utm_factor, "Geographical area", "Cell to utm factor, i.e. utm unit distance between elevation sampling points", Int; kwds...)
    sheet_width_mm = get_kw_or_config_value(:sheet_width_mm ,"Printer consistent capability", "Printable width mm", Int; kwds...)
    sheet_height_mm = get_kw_or_config_value(:sheet_height_mm ,"Printer consistent capability", "Printable height mm", Int; kwds...)
    density_pt_m⁻¹ = get_kw_or_config_value(:density_pt_m⁻¹ ,"Geographical area", "Output density, i.e. 'cells' / 'dots' / 'points' or 'pixels' per paper meter", Int; kwds...)
    pth = get_kw_or_config_value(:pth, "File folder", "Top folders path under homedir()", String; kwds...)
    # This value is for checking if density_pt_m⁻¹ is higher than is printable
    density_limit_pt_inch⁻¹ = get_kw_or_config_value(:density_limit_pt_inch⁻¹ ,"Printer consistent capability", "Stated density limit, dots per inch", Int; kwds...)
    # Some value checks
    @assert density_pt_m⁻¹ <= density_limit_pt_inch⁻¹ / 0.0254 "Printing density exceeds capability, $(Int(round(density_limit_pt_inch⁻¹ / 0.0254)))m⁻¹."
    if isnothing(match(r"(b|B)itmap(M|m)aps", pth))
        throw(ArgumentError("Any BmParition path must match regex r\"(b|B)itmap(M|m)aps\". Current path is $pth"))
    end
    # Either nrc (number of rows and columns) or sheet_indices can be specified.
    # If nrc is specified, this overrules.
    if (:sheet_indices ∉ keys(kwds)) || (:sheet_indices ∈ keys(kwds) && :nrc ∈ keys(kwds))
        nrc = get_kw_or_config_value(:nrc ,"Geographical area", "Output paper sheets (rows columns)", Tuple{Int, Int}; kwds...)
        smb = SheetMatrixBuilder(southwest_corner, nrc, cell_to_utm_factor, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, pth)
    elseif :sheet_indices ∈ keys(kwds) && :nrc ∉ keys(kwds)
        sheet_indices = kwds[:sheet_indices]
        smb = SheetMatrixBuilder(southwest_corner, sheet_indices, cell_to_utm_factor, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, pth)
    else
        throw("Surprise! Why?")
    end
    # Often, the user will be interested in inspecting
    # derived properties rather than the basic ones.
    show_augmented(smb)
    # Nice to know.
    if geo_area(smb) >= 16e6
        @info "Since geographical area per sheet is  > 16km², 'høydedata.no' will not provide a single elevation data file per sheet. The pipeline will make a consolidated single file."
    end
    smb
end


function process_job(smb, complete_sheets_first)
        operations_order = [establish_folder, 
        unzip_tif, 
        consolidate_elevation_data, 
        water_overlay, 
        topo_relief, 
        contour_lines_overlay, 
        grid_overlay, 
        ridge_overlay,
        summit_markers,
        join_layers]
    # 
    if complete_sheets_first
        for sb in smb
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
    # TODO: We did not consider cross-sheet summit prominence yet.
    #       This can not be done before all sheets are processed. 
    #       But the interaction depends on an initial step without
    #       sheet interaction.
    #       We don't want to give up the current flexibility (i.e. 
    #       processing all sheets at once, or finish one first).
    #       Hence, a good approach might be to require re-running
    #       the pipeline after fully finishing. 
    #       Interaction would only be triggered if the required 
    #       files for interaction have been generated.
    #       Since interaction is a SheetMatrixBuilder level operation,
    #       which may take time and require more memory, the prerequisites 
    #       for it should be checked AFTER join_layers, i.e. here.
    true
end
function call_func(fn, sb)
    @debug "Calling `$fn`. $(full_folder_path(sb))"
    ok_res = fn(sb)
    if ! ok_res
        @warn "Could not finish $fn($sb) with success. Exiting"
        return false
    else
        @debug "Finished"
    end
    true
end