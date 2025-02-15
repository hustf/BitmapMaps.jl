# This pipeline is the most central part of the package.
"""
    run_bitmapmap_pipeline(; complete_sheets_first = true, skip_summits = false, kwds...)
    ---> SheetMatrixBuilder

The job is defined in file BitmapMaps.ini, in user's home directory. Call `less(homedir() * "/BitmapMaps.ini")` to inspect.

# Arguments

You can overrule parameters from the .ini file with identical keywords, but changing the .ini file
is recommended. 

- `example_keyword` (would appear as `:example_keyword`) in the file.
- `complete_sheets_first`=true   'true' means one sheet is fully processed, then the next sheet.
                                 'false' means that each operation is finished for all sheets before the next operation.
- `skip_summits`=false           'true' skips the prominence (re-)calculation, which can save time. A few iteration is 
                                 needed for correct boundary conditions between sheets.

# Example
```
julia> run_bitmapmap_pipeline(;nrc = (2, 2));
[ Info: Pipeline running. Set ENV["JULIA_DEBUG"] = "BitmapMaps" for detailed progress.
Bitmapmap builder based on .ini file and keywords
SheetMatrixBuilder((35425, 6920995), # southwest_corner
      CartesianIndices((1:2, 1:2)), # sheet_indices
                                 3, # cell_to_utm_factor
                               191, # sheet_width_mm
                               275, # sheet_height_mm
                             11811, # density_pt_m⁻¹
              "bitmapmaps/default") # pth
        [easting, northing] derived properties:
          Bounding Box (BB) SE-NW            = (35425 6920995)-(48955 6940483)
          Northeast internal corner          = (48952, 6940483) - most northeastern sample point
          Geo centre                         = (42190.0, 6.930739e6)
          Grid centre single                 = (42190, 6930740)
        Derived properties:
          Geographical (width, height) [km]  = (13.5, 19.5)
          Geographical area [km²]            = 264
                    Per sheet = 65.92 km²   (Single file export limit: 16 km²)
          Sheets total (width, height) [cm]  = (38.2, 55.0)
                    Per sheet [mm] (w, h) = (190.9, 275.0)
          Map scale                          = 1 : 35433 = 1 : (cell_to_utm_factor * density_pt_m⁻¹)
        BBs of sheets as Well Known Text (paste in wktmap.com or nvdb-vegdata.github.io/nvdb-visrute/STM ):
          MULTIPOLYGON (
                   ((35425 6920995, 42190 6920995, 42190 6930739, 35425 6930739, 35425 6920995)),
                   ((35425 6930739, 42190 6930739, 42190 6940483, 35425 6940483, 35425 6930739)),
                   ((42190 6920995, 48955 6920995, 48955 6930739, 42190 6930739, 42190 6920995)),
                   ((42190 6930739, 48955 6930739, 48955 6940483, 42190 6940483, 42190 6930739)))
[ Info: No .tif files in C:\\Users\\f\bitmapmaps/default\\1 1  35425 6920995  42190 6930739 to consolidate. Exiting.
[ Info: Could not make Consolidated.tif for sheet  with folder path bitmapmaps/default\1 1  35425 6920995  42190 6930739. Download and unzip .tif files? Exiting.
┌ Warning: Could not finish consolidate_elevation_data(SheetBuilder((0, 3248), (1:3248, 1:2255), (35425, 6930739)@(1, 1), 1, 11811, "bitmapmaps/default\\1 1  35425 6920995  42190 6930739")
│ ) with success. Exiting
```
"""
function run_bitmapmap_pipeline(; complete_sheets_first = true, skip_summits = false, kwds...)
    if get(ENV, "JULIA_DEBUG", "") !== "BitmapMaps"
        @info "Pipeline running. Set ENV[\"JULIA_DEBUG\"] = \"BitmapMaps\" for detailed progress."
    end
    smb = define_builder(; kwds...)
    ok_res = process_job(smb, complete_sheets_first, skip_summits)
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
    allowed_keywords = [:southwest_corner, :cell_to_utm_factor, :sheet_width_mm, :sheet_height_mm, :density_pt_m⁻¹,
        :pth, :density_limit_pt_inch⁻¹, :sheet_indices, :nrc]
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
    smb
end


function process_job(smb, complete_sheets_first, skip_summits)
    # Pre-process: Make an overview .svg file, referring the to-be-established files:
    make_vector_graphics(smb)
    #
    # Closure for making thumbnails. Note, it would be cleaner to add another field to the
    # SheetBuilder type, but making thumbnails is sort of an add-on functionality, so we 
    # don't revise the type (at this late revision).
    n_governing = max(size(smb)...)
    make_thumbnail(sb) = make_thumbnail_image(full_folder_path(sb), sb.density_pt_m⁻¹, n_governing)
    #
    # Prepare closure for summits and prominence.
    # While almost all operations require no interaction between sheets, correct summit prominence
    # typically is affected by neighbouring sheets. We can't store all data for all sheets simultaneously,
    # so we choose to store the relevant data for this calculation in a sparse meta-graph indexed by utm coordinates.
    # Filename for the graph which is common to all sheets
    ffna_graph = joinpath(full_folder_path(smb), bbox_external_string(smb) * ".z")
    @debug """    Common to all sheets, graph file: ffna_graph = raw"$ffna_graph" """
    # Closure
    f_sides_with_border  = let smb = smb
        (sb) -> sides_with_border(smb, sb)
    end
    # Closure
    summits_and_prominence(sb) = summits_on_sheet(sb, ffna_graph, f_sides_with_border)
    # Define the steps order - each step is called with each sheet builder.
    operations_order = [establish_folder,
        unzip_tif,
        consolidate_elevation_data,
        water_overlay,
        topo_relief,
        contour_lines_overlay,
        grid_overlay,
        ridge_overlay,
        summits_and_prominence,
        join_layers,
        make_thumbnail]
    # Do the sheet steps depth-first or width-first.
    if complete_sheets_first
        for sb in smb
            @info "Sheet $(cartesian_index_string(smb, sb.sheet_number)) of up to $(cartesian_index_string(smb))"
            for (i, fn) in enumerate(operations_order)
                call_func(fn, sb, skip_summits, i) || return false
            end
        end
    else
        for (i, fn) in enumerate(operations_order)
            for sb in smb
                @info "Sheet $(cartesian_index_string(smb, sb.sheet_number)) of up to $(cartesian_index_string(smb))"
                call_func(fn, sb, skip_summits, i) || return false
            end
        end
    end
    #
    # Sheet interaction, or 'regional' step. Requires that `summits_on_sheet` have been run for 
    # all sheets first. 
    #
    @info "Summits regional update (find names and prominence)"
    summits_regional_update(smb, ffna_graph)
    #
    @info "Update vector graphics and joining layers for sheets [1, 1] to $(cartesian_index_string(smb))"
    for sb in smb
        call_func(make_vector_graphics, sb, skip_summits, length(operations_order) + 2) || return false
        call_func(join_layers, sb, skip_summits, length(operations_order) + 3) || return false
    end
    true
end
function call_func(fn, sb, skip_summits, i)
    if !skip_summits || fn !== summits_on_sheet
        @debug "$(lpad(i, 2)) `$fn`. $(full_folder_path(sb))"
        ok_res = fn(sb)
        if ! ok_res
            @warn "Could not finish step $i $fn($sb) with success. Exiting"
            return false
        else
            #@debug "Finished `$fn`"
        end
    end
    true
end