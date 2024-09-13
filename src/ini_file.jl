function _prepare_init_file_configuration(io)
    # Add a comment at top (IniFile.jl has no functions for comments)
    msg = """
        # Configuration file for 'BitmapMaps'.
        # You can modify and save the values here. To start over from 'factory settings':
        # Delete this file. A new file will be created next time values stored in the ini file is accessed.
        #
        # This defines the UTM boundaries of your map, and its division into sheets of paper.
        # Sheets of paper will be numbered from South-West (index [1, 1] or [1]) to North-East (index [end, end] or [end]),
        # and a folder structure will be established.
        #
        # Utm coordinates (easting, northing) increase to the East (x) and North (y).
        # A boundary box is defined by two of its external corners:
        #     (minimum_easting minimum_northing)-(maximum_easting maximum_northing)
        # A simple way to get EU89, UTM-zone 33 coordinates is through norgeskart.no.
        # We project geography on the UTM 33 grid. In 'norgeskart.no' that and other
        # grids can be displayed through hamburger button > fastmerker > utm-rutenett.
        #
        """
    println(io, msg)
    # Note, we're creating several 'IniFile' below, for sequential writing to the same file.
    # The sequence keeps the order of heading consistent as planned. We want all those variables 
    # that can't be overruled last.
    # At reading, they are of course parsed as the same Ini File.
    ini = Inifile()
    # Shorthand
    entry(section, key, val; comm = "") = set_section_key_string_comment(ini, section, key, val, comm)
    # Section which can be overruled by keywords to `run_bitmapmap_pipeline`
    entry("Printer consistent capability", "Printable width mm", "191"; comm = ":sheet_width_mm\r\n  #    Measured 193 mm. Allowing 2 mm  random variation.")
    entry("Printer consistent capability", "Printable height mm", "275"; comm = ":sheet_height_mm\r\n  #    Measured 277 mm. Allowing 2 mm for random variation.")
    entry("Printer consistent capability", "Stated density limit, dots per inch", "600"; comm = ":density_limit_pt_inch⁻¹\r\n  #    As advertised by Brother")
    println(io, ini)
    ini = Inifile()
    entry("Geographical area", "Southwest corner (utm easting northing)", "(35425 6920995)"; comm = ":southwest_corner")
    entry("Geographical area", "Output paper sheets (rows columns)", "(3 4)"; comm = ":nrc")
    entry("Geographical area", "Output density, i.e. 'cells' / 'dots' / 'points' or 'pixels' per paper meter", "11811";
        comm = ":density_pt_m⁻¹\r\n  #    For reference, 300  / inch = 300 / (0.0254 m) = 11811 m⁻¹ \r\n  #    Use lower values to slightly 'zoom in', otherwise use 'cell_to_utm_factor'.")
    entry("Geographical area", "Cell to utm factor, i.e. utm unit distance between elevation sampling points", "3"; comm = ":cell_to_utm_factor\r\n  #    How many 'utm metres' does a 'cell' / 'dot' / 'point' or 'pixel' side represent?")
    println(io, ini)
    ini = Inifile()
    entry("File folder", "Top folders path under homedir()", "bitmapmaps/default"; comm = ":pth \r\n  #    Folder hierarchy will be created under homedir().\r\n  #    For copying in files, see 'copy_relevant_tifs_to_folder'")
    println(io, ini)
    msg = """
        #
        # Sections which can not be overruled by keywords to `run_bitmapmap_pipeline`
        #
        """
    println(io, msg)
    # Sections which can not be overruled by keywords to `run_bitmapmap_pipeline`
    ini = Inifile()
    entry("Elevation contour lines", "Thickness 20m", "1  "; comm = "\r\n  #  Values: 0, 1, 3, 5...")
    entry("Elevation contour lines", "Thickness 100m", "3 "; comm = "\r\n  #  Values: 0, 1, 3, 5...")
    entry("Elevation contour lines", "Thickness 1000m", "5"; comm = "\r\n  #  Values: 0, 1, 3, 5...")
    entry("Elevation contour lines", "Minimum length", "6"; comm = "\r\n  #  Mostly affects forests and built-up")
    entry("UTM grid", "Grid line thickness", "5"; comm = "\r\n  #  5 works well for :density_pt_m⁻¹ = 11811")
    entry("UTM grid", "Grid spacing [m]", "1000")
    entry("UTM grid", "Zone for elevation data", "33"; comm = "\r\n  #  File metadata is ignored. wgs84 datum assumed.")
    entry("Water", "Lake steepness max", "0.075"; comm = "\r\n  #  Values with success in different terrains: 0.075, 0.156, 0.16, 0.2")
    entry("Markers", "Prominence level [m], prominent summit", "100")
    entry("Markers", "Prominence level [m], obscure summit", "50")
    entry("Markers", "Symbol prominent summit", "in_triangle"; comm = "
        # on_square, on_triangle, on_circle
        # in_square, in_triangle, in_circle
        # on_cross, on_xcross, on_hline, on_vline")
    entry("Markers", "Symbol obscure summit", "on_triangle")
    entry("Markers", "Size prominent summit symbol", "31"; comm = "\r\n  #  Length of bounding box size. Must be odd.")
    entry("Markers", "Size obscure summit symbol", "21"; comm = "\r\n  #  Length of bounding box size. Must be odd.")
    entry("Markers", "Minimum stress level", "-14"; comm = "\r\n  #  Summits with even lower stress values are discarded as power lines / artifacts")
    entry("Text", "Font size [pt]", "10.0"; comm = "\r\n  #  Unit is postscript points, used in e.g. Office Word\r\n  # pt = inch / 72 = 0.3528 mm")
    entry("Text", "Flip text side at width fraction", "0.85"; comm = "\r\n  #  Set text to the left of summits when summit is further to the right")
    # To file..
    println(io, ini)
end


"""
    get_config_value(sect::String, key::String)
    get_config_value(sect, key, type::DataType; nothing_if_not_found = false)

Instead of passing long argument lists, we store configuration in a text file.
"""
function get_config_value(sect::String, key::String; nothing_if_not_found = false)
    fna = _get_ini_fnam()
    ini = read(Inifile(), fna)
    if sect ∉ keys(sections(ini))
        msg = """$sect not a section in $fna.
        The existing are: $(keys(sections(ini))).
        If you `BitmapMaps.delete_init_file()`, a new template will be generated when needed.
        """
        throw(ArgumentError(msg))
    end
    if nothing_if_not_found
        IniFile.get(ini, sect, key,  nothing)
    else
        s = IniFile.get(ini, sect, key,  "")
        if s == ""
            throw(ArgumentError("""
                $key not a key with value in section $sect of file $fna.

                Example:
                [user]                          # section
                user_name     = slartibartfast  # key and value
                perceived_age = 5

                Current file:
                $ini
            """))
        end
        s
    end
end
function get_config_value(sect, key, type::DataType; nothing_if_not_found = false)
    st = get_config_value(sect, key; nothing_if_not_found)
    isnothing(st) && return nothing
    tryparse(type, split(st, ' ')[1])
end

function get_config_value(sect, key, ::Type{String}; nothing_if_not_found = false)
    st = strip(get_config_value(sect, key; nothing_if_not_found), ' ')
    stv = strip(split(st, '#')[1])
    isnothing(stv) && return nothing
    String(stv)
end

function get_config_value(sect, key, ::Type{Tuple{Int64, Int64}}; nothing_if_not_found = false)
    st = strip(get_config_value(sect, key; nothing_if_not_found), ' ')
    isnothing(st) && return nothing
    stv = split(split(replace(st, "(" => "", ")" => ""), '#')[1], ' ')
    (tryparse(Int64, stv[1]),     tryparse(Int64, stv[2]))
end


"Get an existing, readable ini file name, create it if necessary"
function _get_ini_fnam()
    fna = _get_fnam_but_dont_create_file()
    if !isfile(fna)
        open(_prepare_init_file_configuration, fna, "w+")
        # Launch an editor
        if Sys.iswindows()
            run(`cmd /c $fna`; wait = false)
        end
        println("Default settings stored in $fna")
    end
    fna
end
_get_fnam_but_dont_create_file() =  joinpath(homedir(), "BitmapMaps.ini")



function set_section_key_string_comment(ini::Inifile, section::T, key::T, val::T, comment::T) where T<:String
    if comment == ""
        set(ini, section, key, val)
    else
        set(ini, section, key, val * " # " * comment)
    end
end



function get_kw_or_config_value(sy::Symbol, sect, key, type; kwds...)
    # Check for a common bug:
    if !isempty(kwds)
        if first(kwds)[1] == :kwds
            throw(ArgumentError("Optional keywords: Use splatting in call: ;kwds..."))
        end
    end
    # Keywords, if present, overrule config file values
    if sy ∈ keys(kwds)
        va = kwds[sy]
        if va isa type
            return va
        else
            throw(ArgumentError("Keyword $sy not of type $type"))
        end
    else
        get_config_value(sect, key, type; nothing_if_not_found = false)
    end
end



"""
    delete_init_file()

You will lose manual changes, like your own printer's actual margins.
Next time you need an ini file, default settings are used to make a new one.
"""
function delete_init_file()
    fna = _get_fnam_but_dont_create_file()
    if isfile(fna)
        rm(fna)
        println("Removed $fna")
    else
        println("$fna Didn't and doesn't exist.")
    end
end
