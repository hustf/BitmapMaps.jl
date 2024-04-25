function establish_folder(p::SheetBuilder)
    fullpth = full_folder_path(p)
    if ! ispath(fullpth)
        mkpath(fullpth)
    end
    ispath(fullpth)
end
full_folder_path(p::SheetBuilder) = joinpath(homedir(), p.pthsh)
full_folder_path(p::SheetMatrixBuilder) = joinpath(homedir(), p.pth)

function parse_folder_name(fofo)
    @show fofo
    locfo = split(fofo, '\\')[end]
    ints = split(locfo)
    r, c, min_x, min_y, max_x, max_y = parse.(Int, ints)
    @assert min_x < max_x
    @assert min_y < max_y
    r, c, min_x, min_y, max_x, max_y
end