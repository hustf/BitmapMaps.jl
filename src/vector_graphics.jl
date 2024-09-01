# Step in pipeline.
# This file copies and modifies an existing .svg file,
# so as to include by reference the bitmap,
# and adds text for the foreground, based on 'Summits.csv'.
# The advantage to this approach over the previous one
# (LuxorLayout) is scriptable, searchable graphics files,
# and we avoid depending on the Cairo / Pango library.

"""
    make_vector_graphics(sb::SheetBuilder)
    make_vector_graphics(fofo, cell_iter, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹)
    --> Bool

"""
function make_vector_graphics(sb::SheetBuilder)
    make_vector_graphics(full_folder_path(sb), sb.cell_iter, width_adjusted_mm(sb), height_adjusted_mm(sb), sb.density_pt_m⁻¹)
end
function make_vector_graphics(fofo, cell_iter, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹)
    #
    # Early exits
    #
    if ! isfile(joinpath(fofo, COMPOSITE_FNAM))
        @debug "    $COMPOSITE_FNAM in $fofo\n           does not exist. Exiting ´make_vector_graphics´"
        return false
    end
    bitmap_time = mtime(joinpath(fofo, COMPOSITE_FNAM))
    ffna_csv = joinpath(fofo, SUMMITS_FNAM)
    if ! isfile(ffna_csv)
        @debug "    $SUMMITS_FNAM in $fofo\n           does not exist. Exiting ´make_vector_graphics´"
        return false
    end
    summits_time = mtime(ffna_csv)
    ffna_svg = replace(joinpath(fofo, COMPOSITE_FNAM), ".png" => ".svg")
    ffna_css = replace(joinpath(fofo, COMPOSITE_FNAM), ".png" => ".css")
    if isfile(ffna_svg) && isfile(ffna_css)
        svg_time = mtime(ffna_svg)
        if svg_time > bitmap_time && svg_time > summits_time
            @debug "    $(splitpath(ffna_svg)[end]) in $fofo\n           already exists and is newer than its input files. Exiting ´make_vector_graphics´"
            return true
        end
    end
    # Fetch font size and left/ right alignment limit from .ini file
    font_pt = get_config_value("Text", "Font size [pt]", Float32)
    limit_left_align = get_config_value("Text", "Flip text side at width fraction", Float32)
    # Do the work
    _make_vector_graphics(ffna_svg, ffna_css, ffna_csv, cell_iter, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, font_pt, limit_left_align)
    true
end

"""
    _make_vector_graphics(ffna_svg, ffna_css, ffna_csv, cell_iter, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, font_pt, limit_left_align)
    --> Image

- Copy template .svg and .css into fofo, under new names.
- Interpolate dimensions and line distance into .svg file, and font size into css file.
- Read summits file and generate text element.

TODO make inclusion of prominence configureable.
"""
function _make_vector_graphics(ffna_svg, ffna_css, ffna_csv, cell_iter, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, font_pt, limit_left_align)
    #
    # Prepare arguments
    #
    fontsize_px = Int(round(font_pt * 0.0003528 * density_pt_m⁻¹))
    lineheight_px = Int(round(1.2 * fontsize_px))
    # We trust the column order, as defined in `write_prominence_to_csv` and `add_names_to_csv`
    summits_data = readdlm(ffna_csv, '\t')[2:end,:]
    if isempty(summits_data)
        @debug "    No summits in file => no text to add."
        return String[]
    end
    vz = Int.(summits_data[:, 1])
    vpr = Int.(round.(summits_data[:, 2]))
    vname = strip.(string.(summits_data[:, 6]))
    v_index = strip.(string.(summits_data[:, 4]))
    vI = Tuple.(split.(map(s -> strip(s,['(', ')']), v_index), ','))
    nx = size(cell_iter, 2)
    vtext_on_left = map(vI) do (_, x)
        tryparse(Float32, x) > (limit_left_align * nx)
    end
    #
    # Do the work
    #
    copy_templates_to_folder(ffna_svg, ffna_css)
    modify_css_font_size(ffna_css, fontsize_px)
    modify_svg_to_sheet_size(ffna_svg, cell_iter, sheet_width_mm, sheet_height_mm)
    modify_svg_text(ffna_svg, vz, vpr, vname, vI, "fjell", lineheight_px, vtext_on_left)
end

function modify_svg_to_sheet_size(ffna_svg, cell_iter, sheet_width_mm, sheet_height_mm)
    ny, nx = size(cell_iter)
    doc = readxml(ffna_svg)
    ns = ["x" => "http://www.w3.org/2000/svg"]
    svg = EzXML.root(doc)
    svg["width"] = "$(round(sheet_width_mm, digits=3))mm"
    svg["height"] = "$(round(sheet_height_mm, digits=3))mm"
    svg["viewBox"] = "1 1 $(nx - 1) $(ny - 1)"
    xp = "x:defs/x:filter/x:feImage"
    feImage = findfirst(xp, svg, ns)
    feImage["width"] = "$(nx - 1)"
    feImage["height"] = "$(ny - 1)"
    commentnode = findfirst("x:defs / x:filter / comment()", svg, ns)
    setnodecontent!(commentnode, """
        Ensure the image fills the rectangle and is positioned correctly.
        The values in this svg file match the size of $COMPOSITE_FNAM with
        (w, h) = ($nx, $ny) pixels.
    """)
    write(ffna_svg, doc)
end

function modify_svg_text(ffna_svg, vz, vpr, vname, vI, class, lineheight_px, vtext_on_left)
    #
    # Get the xml structure into memory
    #
    doc = readxml(ffna_svg)
    ns = ["x" => "http://www.w3.org/2000/svg"]
    svg = EzXML.root(doc)
    # Remove the dummy text element in the template.
    unlink!(findfirst("x:text", svg, ns))
    # Prepare arguments
    dx = 0.2 * lineheight_px
    #
    # Add summit text, one at a time
    #
    for (name, z, prominence, (y, x), text_on_left) in zip(vname, vz, vpr, vI, vtext_on_left)
        if name == ""
            add_single_line_text_element(svg, "$(z)m  ($(prominence))", class,strip(x), strip(y), dx, text_on_left)
        else
            add_two_line_text_element(svg, name, "$(z)m  ($(prominence))", class,strip(x), strip(y), lineheight_px, dx, text_on_left)
        end
    end
    #
    # Save xml
    #
    write(ffna_svg, doc)
    nothing
end
function add_two_line_text_element(parent, text1, text2, class, x, y, lineheight_px, dx, text_on_left)
    # We're misusing indent to reflect xml structure in this function
        el = ElementNode("text")
        link!(el, AttributeNode("class", class))
        link!(el, AttributeNode("x", x))
        link!(el, AttributeNode("y", y))
        link!(el, TextNode("\n    ")) # line break and tabs for xml readability
            tsp1 = ElementNode("tspan")
            link!(tsp1, TextNode(text1))
            if text_on_left
                link!(tsp1, AttributeNode("dx", string(-dx)))
                link!(tsp1, AttributeNode("text-anchor", "end"))
            else
                link!(tsp1, AttributeNode("dx", string(dx)))
            end
            link!(tsp1, TextNode("\n    ")) # line break and tabs for xml readability
        link!(el, tsp1)
            tsp2 = ElementNode("tspan")
            link!(tsp2, TextNode(text2))
            link!(tsp2, AttributeNode("x", x))
            link!(tsp2, AttributeNode("y", y))
            link!(tsp2, AttributeNode("dy", "$lineheight_px"))
            if text_on_left
7                link!(tsp2, AttributeNode("dx", string(-dx)))
                link!(tsp2, AttributeNode("text-anchor", "end"))
            else
                link!(tsp2, AttributeNode("dx", string(dx)))
            end
            link!(tsp2, TextNode("\n    ")) # line break and tabs for xml readability
        link!(el, tsp2)
    link!(parent, el)
    # Add a line break after the </text>, for readability.
    link!(parent, TextNode("\n"))
end
function add_single_line_text_element(parent, text, class, x, y, dx, text_on_left)
    el = ElementNode("text")
    link!(el, TextNode(text))
    link!(el, AttributeNode("class", class))
    link!(el, AttributeNode("x", x))
    link!(el, AttributeNode("y", y))
    if text_on_left
        link!(el, AttributeNode("dx", string(-dx)))
        link!(el, AttributeNode("text-anchor", "end"))
    else
        link!(el, AttributeNode("dx", string(dx)))
    end
    link!(parent, el)
    # Add a line break after the </text>, for readability.
    link!(parent, TextNode("\n"))
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