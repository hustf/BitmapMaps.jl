# Step in pipeline.
# Called prior to the sheets pipeline with a SheetMatrixBuilder to make a navigable .svg file.
# This file copies and modifies an existing .svg file,
# so as to include by reference the bitmap,
# and adds text for the foreground, based on 'Summits.csv'.
# The advantage to this approach over the previous one
# (LuxorLayout) is scriptable, searchable graphics files,
# and we avoid depending on the Cairo / Pango library.
#
# We make vector graphics files because it works better with searchable text.
# A bitmap is referred (included in) each vector graphics.

"""
    make_vector_graphics(smb::SheetMatrixBuilder)
    make_vector_graphics(sb::SheetBuilder; revise_names = false))
    make_vector_graphics(fofo, cell_iter, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹; revise_names = false))
    --> Bool

Create an svg file.

- When called with a SheetMatrixBuilder, makes an .svg mosaic from the sheets.
The mosaic tiles show the topographic relief without overlays, but links to the composite images.
- When called with a SheetBuilder, makes an .svg with the composite bitmap overlain with vector graphics text.
"""
function make_vector_graphics(smb::SheetMatrixBuilder)
    fna_svg = bbox_external_string(smb) * ".svg"
    ffna_svg = joinpath(full_folder_path(smb), fna_svg)
    @info "Mosaic of linked sheets thumbnails in $ffna_svg"
    # Early exit
    if isfile(ffna_svg)
        @debug "    $ffna_svg in $(full_folder_path(smb))\n           already exists. Exiting `make_vector_graphics`"
        return true
    end
    # Do the work
    make_vector_mosaic(smb, ffna_svg)
end
function make_vector_graphics(sb::SheetBuilder; revise_names = false)
    make_vector_graphics(full_folder_path(sb), sb.cell_iter, width_adjusted_mm(sb), height_adjusted_mm(sb), sb.density_pt_m⁻¹; revise_names)
end
function make_vector_graphics(fofo, cell_iter, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹; revise_names = false)
    #
    # Early exits
    #
    # Note, we might continue just fine without Composite.png. That behaviour might be hard to reason about, though.
    if ! isfile(joinpath(fofo, COMPOSITE_FNAM))
       @debug "    $COMPOSITE_FNAM in $fofo\n           does not exist. Exiting ´make_vector_graphics´"
       return false
    end
    bitmap_time = mtime(joinpath(fofo, COMPOSITE_FNAM))
    ffna_csv_summits = joinpath(fofo, SUMMITS_FNAM)
    if ! isfile(ffna_csv_summits)
        @debug "    $SUMMITS_FNAM in $fofo\n           does not exist. Exiting ´make_vector_graphics´"
        return false
    end
    #
    ffna_csv_lakes = splitext(joinpath(fofo, WATER_FNAM))[1] * ".csv"
    if ! isfile(ffna_csv_lakes)
        @debug "    $(splitext(WATER_FNAM)[1] * ".csv") in $fofo\n           does not exist. Exiting ´make_vector_graphics´"
        return false
    end
    #
    summits_time = mtime(ffna_csv_summits)
    lakes_time = mtime(ffna_csv_lakes)
    ffna_svg = replace(joinpath(fofo, COMPOSITE_FNAM), ".png" => ".svg")
    ffna_css = replace(joinpath(fofo, COMPOSITE_FNAM), ".png" => ".css")
    if isfile(ffna_svg) && isfile(ffna_css)
        svg_time = mtime(ffna_svg)
        if svg_time > bitmap_time && svg_time > summits_time && svg_time > lakes_time && ! revise_names
            @debug "    $(splitpath(ffna_svg)[end]) in $fofo\n           already exists and is newer than its input files. Exiting ´make_vector_graphics´"
            return true
        end
    end
    #
    # Fetch font size and left/ right alignment limit from .ini file
    #
    font_pt = get_config_value("Text", "Font size [pt]", Float32)
    limit_left_align = get_config_value("Text", "Flip text side at width fraction", Float32)
    # Do the work
    _make_vector_graphics(ffna_svg, ffna_css, ffna_csv_summits, ffna_csv_lakes, cell_iter, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, font_pt, limit_left_align)
    true
end

"""
    _make_vector_graphics(ffna_svg, ffna_css, ffna_csv_summits, ffna_csv_lakes, cell_iter, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, font_pt, limit_left_align)
    --> Int, creates ffna_svg

- Copy template .svg and .css into fofo, under new names.
- Interpolate dimensions and line distance into .svg file, and font size into css file.
- Read summits file and generate text element.
- Read lakes file and generate text elements.

TODO make inclusion of prominence configureable.
"""
function _make_vector_graphics(ffna_svg, ffna_css, ffna_csv_summits, ffna_csv_lakes, cell_iter, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, font_pt, limit_left_align)
    # Early exit
    if isempty(readdlm(ffna_csv_summits, '\t')[2:end,:]) && isfile(ffna_svg)
        @debug "    No summits in csv file => no text to add."
        return
    end
    #
    # Prepare arguments
    #
    fontsize_px = Int(round(font_pt * 0.0003528 * density_pt_m⁻¹))
    lineheight_px = Int(round(1.2 * fontsize_px))
    nx = size(cell_iter, 2)
    max_x_left_align = limit_left_align * nx
    #
    # Do the work
    #
    copy_templates_to_folder(ffna_svg, ffna_css)
    modify_css_font_size(ffna_css, fontsize_px)
    modify_svg_to_sheet_size(ffna_svg, ffna_css, cell_iter, sheet_width_mm, sheet_height_mm)
    modify_svg_text(ffna_svg, ffna_csv_summits, ffna_csv_lakes, lineheight_px, max_x_left_align)
end

function modify_svg_to_sheet_size(ffna_svg, ffna_css, cell_iter, sheet_width_mm, sheet_height_mm)
    ny, nx = size(cell_iter)
    doc = readxml(ffna_svg)
    ns = ["x" => "http://www.w3.org/2000/svg"]
    # Update style sheet ref.
    stylenode = firstnode(doc.node)
    stylecontent = "type=\"text/css\" href=\"$(splitpath(ffna_css)[end])\""
    setnodecontent!(stylenode, stylecontent)
    # Update svg contents
    svg = EzXML.root(doc)
    svg["width"] = "$(round(sheet_width_mm, digits=3))mm"
    svg["height"] = "$(round(sheet_height_mm, digits=3))mm"
    svg["viewBox"] = "1 1 $(nx - 1) $(ny - 1)"
    xp = "x:defs/x:filter/x:feImage"
    feImage = findfirst(xp, svg, ns)
    feImage["width"] = "$(nx - 1)"
    feImage["height"] = "$(ny - 1)"
    commentnode = findfirst("x:defs / comment()", svg, ns)
    setnodecontent!(commentnode, """
        Ensure the image fills the rectangle and is positioned correctly.
        The values in this svg file match the size of $COMPOSITE_FNAM with
        (w, h) = ($nx, $ny) pixels.
    """)
    # This error occured two times at start-up:
    #
    # ┌ Debug: Pre-process: Make C:\Users\f\BitmapMaps/Stetind\(564923 7560195)-(567177 7563442).svg
    # └ @ BitmapMaps C:\Users\f\.julia\packages\BitmapMaps\jlaPp\src\vector_graphics.jl:339
    # ERROR: XMLError: Permission denied: C:\Users\f\BitmapMaps/Stetind\(564923 7560195)-(567177 7563442).svg from Input/Output stack (code: 1501, line: 0)
    # Stacktrace:
    #   [1] throw_xml_error()
    #     @ EzXML C:\Users\f\.julia\packages\EzXML\qbIRq\src\error.jl:87
    #   [2] macro expansion
    #     @ C:\Users\f\.julia\packages\EzXML\qbIRq\src\error.jl:52 [inlined]
    #   [3] write(filename::String, doc::EzXML.Document)
    #     @ EzXML C:\Users\f\.julia\packages\EzXML\qbIRq\src\document.jl:247
    #   [4] modify_svg_to_sheet_size(ffna_svg::String, ffna_css::String, cell_iter::CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}, sheet_width_mm::Int64, sheet_height_mm::Int64)
    #     @ BitmapMaps C:\Users\f\.julia\packages\BitmapMaps\jlaPp\src\vector_graphics.jl:139
    #
    #
    # Since we are not able to re-create, we suspect it has something to do with compilation of EzXML,
    # and desperately insert an arbitrary wait here:
    sleep(1)
    # This is where an error might trigger:
    write(ffna_svg, doc)
end

function modify_svg_text(ffna_svg, ffna_csv_summits, ffna_csv_lakes, lineheight_px, max_x_left_align)
    #
    # Get the xml structure into memory
    #
    doc = readxml(ffna_svg)
    ns = ["x" => "http://www.w3.org/2000/svg"]
    svg = EzXML.root(doc)
    # Remove the dummy text element in the template.
    unlink!(findfirst("x:text", svg, ns))
    # Add summit and lakes text
    add_summit_text(svg, ffna_csv_summits, lineheight_px, max_x_left_align)
    add_lake_text(svg, ffna_csv_lakes)
    #
    # Save xml
    #
    write(ffna_svg, doc)
    nothing
end
function add_summit_text(svg, ffna_csv_summits, lineheight_px, max_x_left_align)
    # Prepare
    class = "fjell"
    subclass = "fjell_alt"
    dx = 0.2 * lineheight_px
    # Read the text contents
    # We trust the column order, as defined in `write_prominence_to_csv` and `add_names_to_csv`
    summits_data = readdlm(ffna_csv_summits, '\t')[2:end,:]
    if isempty(summits_data)
        @debug "    No summits in file $(joinpath(splitpath(ffna_csv_summits)[end - 1: end])) \n\t\t => no text to add."
        return
    end
    if size(summits_data, 2) < 6
        @debug "    No summit names in file $(joinpath(splitpath(ffna_csv_summits)[end - 1: end])) \n\t\t => no text to add."
        return
    end
    vz = Int.(summits_data[:, 1])
    vpr = Int.(round.(summits_data[:, 2]))
    vname = strip.(string.(summits_data[:, 6]))
    v_index = strip.(string.(summits_data[:, 4]))
    vI = Tuple.(split.(map(s -> strip(s,['(', ')']), v_index), ','))
    vtext_on_left = map(vI) do (_, x)
        tryparse(Float32, x) > max_x_left_align
    end
    #
    # Add summit text, one at a time
    #
    for (name, z, prominence, (y, x), text_on_left) in zip(vname, vz, vpr, vI, vtext_on_left)
        text_anchor = text_on_left ? "end" : "start"
        if name == ""
            inject_text!(svg, "$(z)m  ($(prominence)m)", strip(x), strip(y), dx, text_anchor; class, subclass)
        else
            inject_two_line_text!(svg, name, "$(z)m  ($(prominence)m)", strip(x), strip(y), lineheight_px, dx, text_anchor; class, subclass)
        end
    end
end

function add_lake_text(svg, ffna_csv_lakes)
    # Read the text contents
    # We trust the column order, not checking headers
    lakes_data = readdlm(ffna_csv_lakes, '\t')[2:end,:]
    if isempty(lakes_data)
        @debug "    No lakes in file => no text to add."
        return
    end
    vz = Int.(lakes_data[:, 1])
    v_index = strip.(string.(lakes_data[:, 4]))
    vI = Tuple.(split.(map(s -> strip(s,['(', ')']), v_index), ','))
    #
    # Add lake text, one at a time
    #
    for (z, (y, x)) in zip(vz, vI)
        if z > 1
            inject_text!(svg, "$(z)m", strip(x), strip(y), 0, "middle"; class = "sjø")
        end
    end
end


function inject_two_line_text!(parent, text1, text2, x, y, lineheight_px, dx, text_anchor; class = nothing, subclass = nothing, tab_lev = 0)
    el = ElementNode("text")
    isnothing(class) || link!(el, AttributeNode("class", class))
    link!(el, AttributeNode("x", x))
    link!(el, AttributeNode("y", y))
    inject_tspan!(el, text1; dx, text_anchor, subclass)
    inject_tspan!(el, text2; dx, text_anchor, x, y, lineheight_px, subclass)
    link!(parent, el)
    link!(parent, TextNode('\n' * repeat("    ", tab_lev)))  # line break and indent for xml readability
    nothing
end

function inject_text!(parent, text, x, y, dx, text_anchor; class = nothing, subclass = nothing, tab_lev = 0)
    el = ElementNode("text")
    isnothing(class) || link!(el, AttributeNode("class", class))
    link!(el, AttributeNode("x", x))
    link!(el, AttributeNode("y", y))
    add_alignment!(el, dx, text_anchor)
    #link!(el, TextNode('\n' * repeat("    ", tab_lev + 1))) # line break and indent for xml readability
    #link!(el, TextNode(text))
    link!(el, TextNode('\n' * repeat("    ", tab_lev + 1))) # line break and indent for xml readability
    # Separate enclosed paranthesis, if any
    rgx = r"\((.+?)\)"
    m = match(rgx, text)
    if isnothing(m)
        link!(el, TextNode(text))
    else
        # Add the text outside of paranthesis to this element, tsp
        link!(el, TextNode(replace(text, rgx => "")))
        # Inject another tspan with the paranthesis content into this one
        inject_tspan!(el, first(m.captures); class = subclass, tab_lev = tab_lev + 1)
    end
    #
    link!(el, TextNode('\n' * repeat("    ", tab_lev))) # line break and indent for xml readability
    link!(parent, el)
    link!(parent, TextNode('\n' * repeat("    ", tab_lev)))  # line break and indent for xml readability
    nothing
end

"""
    inject_tspan!(parent, text; dx = nothing, text_anchor = nothing,
        x = nothing, y = nothing, lineheight_px = nothing,
        class = nothing, subclass = nothing, tab_lev = 1)

# Arguments

- 'class' applies to the element itself
- 'subclass' applies to the recursive call when enclosed paranthesises are part of 'text'.
- 'text' is the text string. If it contains enclosed paranthesis, those are removed and
   a tspan is injected within the tspan element through a recursive call.
"""
function inject_tspan!(parent, text; dx = nothing, text_anchor = nothing,
    x = nothing, y = nothing, lineheight_px = nothing,
    class = nothing, subclass = nothing, tab_lev = 1)
    #
    link!(parent, TextNode('\n' * repeat("    ", tab_lev))) # line break and indent for xml readability
    tsp = ElementNode("tspan")
    isnothing(x) || link!(tsp, AttributeNode("x", x))
    isnothing(y) || link!(tsp, AttributeNode("y", y))
    isnothing(lineheight_px) || link!(tsp, AttributeNode("dy", "$lineheight_px"))
    isnothing(class) || link!(tsp, AttributeNode("class", class))
    if ! isnothing(dx)
        if ! isnothing(text_anchor)
            add_alignment!(tsp, dx, text_anchor)
        end
    end
    link!(tsp, TextNode('\n' * repeat("    ", tab_lev + 1))) # line break and indent for xml readability
    # Separate enclosed paranthesis, if any
    rgx = r"\((.+?)\)"
    m = match(rgx, text)
    if isnothing(m)
        link!(tsp, TextNode(text))
    else
        # Add the text outside of paranthesis to this element, tsp
        link!(tsp, TextNode(replace(text, rgx => "")))
        # Inject another tspan with the paranthesis content into this one, recursively.
        inject_tspan!(tsp, first(m.captures); class = subclass, tab_lev = tab_lev + 1)
    end
    link!(tsp, TextNode('\n' * repeat("    ", tab_lev))) # line break and indent for xml readability
    link!(parent, tsp)
end

function add_alignment!(el, dx, text_anchor::String)
    if text_anchor == "end" || text_anchor == "middle"
        link!(el, AttributeNode("dx", string(-dx)))
        link!(el, AttributeNode("text-anchor", text_anchor))
    elseif text_anchor == "start"
        link!(el, AttributeNode("dx", string(dx)))
    else
        throw(ArgumentError("text_anchor: $text_anchor"))
    end
end


function copy_templates_to_folder(ffna_svg, ffna_css)
    ffna_template_svg = joinpath(@__DIR__, "..", "resource", "template.svg")
    ffna_template_css = joinpath(@__DIR__, "..", "resource", "template.css")
    isfile(ffna_template_svg) || throw(ErrorException("Could not find file $ffna_template_svg"))
    isfile(ffna_template_css) || throw(ErrorException("Could not find file $ffna_template_css"))
    # We will always overwrite existing files
    cp(ffna_template_svg, ffna_svg; force = true)
    cp(ffna_template_css, ffna_css; force = true)
end

function modify_css_font_size(ffna_css, fontsize_px)
    vs = readlines(ffna_css, keep = true)
    vs[2] = replace(vs[2], "fontsize: 78px;" => "fontsize: $(fontsize_px)px;")
    open(ffna_css, "w") do io
        foreach(vs) do s
            write(io, s)
        end
    end
end

"""
    make_vector_mosaic(smb, ffna_svg)

Called via make_vector_graphics. We refer files not knowing if they're made yet.
"""
function make_vector_mosaic(smb, ffna_svg)
    @debug "Pre-process: Make $(ffna_svg)"
    # The mosaic svg currently don't use features 
    # from the .css style sheet, but
    # that may change in future. So make both .svg and .css
    ffna_css = splitext(ffna_svg)[1] * ".css"
    copy_templates_to_folder(ffna_svg, ffna_css)
    println("jj--------------------------------------------------")
    modify_svg_to_sheet_size(ffna_svg, ffna_css, smb[1].cell_iter, smb.sheet_width_mm, smb.sheet_height_mm)
    println("kk--------------------------------------------------")

    make_reference_mosaic(ffna_svg, smb)
end
function make_reference_mosaic(ffna_svg, smb)
    doc = readxml(ffna_svg)
    ns = ["x" => "http://www.w3.org/2000/svg"]
    svg = EzXML.root(doc)
    # Remove the dummy elements in the template.
    unlink!(findfirst("x:text", svg, ns))
    unlink!(findfirst("x:defs / x:filter", svg, ns))
    unlink!(findfirst("x:rect", svg, ns))
    # Prepare params which are identical to every sheet tile
    ny, nx = size(smb[1].cell_iter)
    n_governing = max(size(smb)...)
    tile_width = Int(floor((nx - 1) / n_governing))
    tile_height = Int(floor((ny - 1) / n_governing))
    # Add sheet tiles
    for sb in smb
        # Windows has no problems interpreting '/' as '\'. Here, we're making urls,
        # so we use '/'.
        imgpath = replace(joinpath(splitpath(sb.pthsh)[end], THUMBNAIL_FNAM), "\\" => "/")
        urlpath = replace(joinpath(splitpath(sb.pthsh)[end], replace(COMPOSITE_FNAM, ".png" => ".svg")), "\\" => "/")
        r, c = row_col_of_sheet(smb, sb.sheet_number)
        x, y = Int.(round.(sb.pixel_origin_ref_to_bitmapmap ./ n_governing))
        add_sheet_tile!(svg, x, y, tile_width, tile_height, r, c, imgpath, urlpath, ns)
    end
    write(ffna_svg, doc)
    nothing
end
function add_sheet_tile!(svg, x, y, tile_width, tile_height, r, c, imgpath, urlpath, ns)
    defs = findfirst("x:defs", svg, ns)
    add_filter_def!(defs, x, y, tile_width, tile_height, r, c, imgpath, ns)
    add_linked_tile!(svg, x, y, tile_width, tile_height, r, c, urlpath, ns)
    svg
end

function add_filter_def!(parent, x, y, tile_width, tile_height, r, c, imgpath, ns)
    el = ElementNode("filter")
    link!(el, AttributeNode("id", "f_$(r)_$(c)"))
    link!(el, TextNode("\n\t\t")) # Line break for readability.
    elc1 = ElementNode("feImage")
    link!(elc1, AttributeNode("href", imgpath))
    link!(elc1, AttributeNode("x", "$(x + 0.5)"))
    link!(elc1, AttributeNode("y", "$(y + 0.5)"))
    link!(elc1, AttributeNode("width", "$(tile_width)"))
    link!(elc1, AttributeNode("height", "$(tile_height)"))
    link!(el, elc1)
    link!(el, TextNode("\n\t\t")) # Line break for readability.
    elc2 = ElementNode("feMerge")
    elc2c = ElementNode("feMergeNode")
    link!(elc2c, AttributeNode("in", "img"))
    link!(elc2, elc2c)
    link!(el, elc2)
    link!(el, TextNode("\n\t")) # Line break for readability.
    link!(parent, TextNode("\n\t")) # Line break for readability.
    link!(parent, el)
    link!(parent, TextNode("\n\t")) # Line break for readability.
    parent
end
function add_linked_tile!(parent, x, y, tile_width, tile_height, r, c, urlpath, ns)
    el = ElementNode("a")
    link!(el, AttributeNode("href", urlpath))
    link!(el, TextNode("\n\t\t\t")) # Line break for readability.
    elc = ElementNode("rect")
    link!(elc, AttributeNode("x", "$(x)"))
    link!(elc, AttributeNode("y", "$(y)"))
    link!(elc, AttributeNode("width", "$(tile_width)"))
    link!(elc, AttributeNode("height", "$(tile_height)"))
    fref = "filter:url(#f_$(r)_$(c));"
    link!(elc, AttributeNode("style", fref))
    link!(el, elc)
    link!(parent, TextNode("\n")) # Line break for readability.
    link!(parent, el) # Line break for readability.
    link!(parent, TextNode("\n")) # Line break for readability.
end