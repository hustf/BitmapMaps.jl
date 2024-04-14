
function _prepare_init_file_configuration(io)
    # Add a comment at top (IniFile.jl has no functions for comments)
    msg = """
        # Configuration file for 'BitmapMaps'.
        # You can modify and save the values here. To start over from 'factory settings':
        # Delete this file. A new file will be created next time values stored in the ini file is accessed.
        #
        # This defines the UTM borders of your map, and it's mapping to sheets of paper.
        # Sheets of paper will be numbered from North-West (1, 1) to South-East (n, m), 
        # and a folder structure will be established if necessary.
        #
        # Utm coordinates (easting, northing) increase to the east (x) and north (y). A boundary box is
        # defined by its corners: (minimum_easting minimum_northing)-(maximum_easting maximum_northing)
        # A simple way to get EU89, UTM-zone 33 coordinates is through norgeskart.no.
        # We also project geography on the UTM 33 grid. In 'norgeskart.no' that and other 
        # grids can be displayed through hamburger > fastmerker > utm-rutenett.
        #
        """
    println(io, msg)
    #
    ini = Inifile()
    # Shorthand
    entry(section, key, val; comm = "") = set_section_key_string_comment(ini, section, key, val, comm)
    # Entries in arbitrary order from file. Will be ordered by section.
    entry("Geographical position", "Southwest corner (easting northing)", "(4873 6909048)"; comm = ":southwest_corner")
    # Printed and measured 'resource/printer_test_gutter' with 'fill page' scaling settings
    entry("Printer consistent capability", "Printable width mm", "191"; comm = ":pwi Measured 193, 2 mm for random vari.")
    entry("Printer consistent capability", "Printable height mm", "275"; comm = ":phe Measured 277, 2 mm for random vari.")
    entry("Printer consistent capability", "Stated density limit, dots per inch", "600"; comm = ":pdensmax As advertised by Brother")
    entry("Printing pixel density", "Selected density, dots per inch", "300"; comm = ":pdens Based on assumed viewing distance > 20 cm")
    entry("Number of printable sheets", "(rows columns)", "(3 4)"; comm = ":nrc")
    entry("Pixel to utm factor", "Pixel distance between elevation sampling points", "3"; comm = ":pix_to_utm_factor One pixel / point distance equals 3 metres easting or northing")
    # To file..
    println(io, ini)
end


"""
    get_config_value(sect::String, key::String)
    get_config_value(sect, key, type::DataType; nothing_if_not_found = false)

Instead of passing long argument lists, we store configuration in a text file.
"""
function get_config_value(sect::String, key::String; nothing_if_not_found = false)
    fnam = _get_ini_fnam()
    ini = read(Inifile(), fnam)
    if sect ∉ keys(sections(ini))
        msg = """$sect not a section in $fnam.
        The existing are: $(keys(sections(ini))).
        If you delete the .ini file above, a new template will be generated.
        """
        throw(ArgumentError(msg))
    end
    if nothing_if_not_found
        get(ini, sect, key,  nothing)
    else
        s = get(ini, sect, key,  "")
        if s == ""
            throw(ArgumentError("""
                $key not a key with value in section $sect of file $fnam.

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

