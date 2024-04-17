function establish_folder(p::SheetPartition)
    fullpth = joinpath(homedir(), p.pthsh)
    if ! ispath(fullpth)
        mkpath(fullpth)
    end
    ispath(fullpth)
end