# Step in pipeline.
"""
establish_folder(p::SheetBuilder) ---> Bool

This also saves a Julia parseable file with the sheet's definitition.
"""
function establish_folder(p::SheetBuilder)
    fullpth = full_folder_path(p)
    if ! ispath(fullpth)
        mkpath(fullpth)
    end
    #
    filename = joinpath(fullpth, PARSEABLE_FNAM)
    #
    open(filename, "w") do io
        show(io, MIME"text/plain"(), p)
    end
    ispath(fullpth)
end

full_folder_path(p::SheetBuilder) = joinpath(homedir(), p.pthsh)
full_folder_path(p::SheetMatrixBuilder) = joinpath(homedir(), p.pth)

function parse_folder_name(fofo)
    locfo = split(fofo, '\\')[end]
    if ! startswith(locfo, r"[0-9]")
        throw(ArgumentError("Folder name \"$locfo\" does not match the integer naming scheme: \"r c min_x min_y max_x max_y\"" ))
    end
    ints = split(locfo)
    r, c, min_x, min_y, max_x, max_y = parse.(Int, ints)
    @assert min_x < max_x
    @assert min_y < max_y
    r, c, min_x, min_y, max_x, max_y
end