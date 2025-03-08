# Step in pipeline. Create directories for one sheet.
# Also related utility functions.

"""
    establish_folder(sb::SheetBuilder)
    ---> Bool

This also saves a Julia parseable file with the sheet's definitition.
"""
function establish_folder(sb::SheetBuilder)
    fullpth = full_folder_path(sb)
    if ! ispath(fullpth)
        mkpath(fullpth)
    end
    #
    filename = joinpath(fullpth, PARSEABLE_FNAM)
    #
    open(filename, "w") do io
        show(io, MIME"text/plain"(), sb)
    end
    ispath(fullpth)
end

full_folder_path(sb::SheetBuilder) = joinpath(homedir(), sb.pthsh)
full_folder_path(smb::SheetMatrixBuilder) = joinpath(homedir(), smb.pth)

"""
    parse_folder_name(fofo)
    ---> NTuple{6, Int64}
"""
function parse_folder_name(fofo)
    locfo = split(fofo, '\\')[end]
    if ! startswith(locfo, r"[0-9]")
        throw(ArgumentError("Folder name \"$locfo\" does not match the integer naming scheme: \"r-c_minx-miny_maxx-maxy\"" ))
    end
    ints = split(locfo, r"[_-]", keepempty = false)
    r, c, min_x, min_y, max_x, max_y = parse.(Int, ints)
    @assert min_x < max_x
    @assert min_y < max_y
    r, c, min_x, min_y, max_x, max_y
end