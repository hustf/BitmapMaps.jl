
"""
    unzip_tif(sb::SheetBuilder)
    unzip_tif(fofo)
    ---> Bool

Unzip the .tif files in downloaded zip. Typically, unzipped files are
largish, but easily compressable. See `tif_full_filenames_buried_in_folder`
for candidates for deletion. Just keep the zip file and CONSOLIDATED_FNAM.
"""
unzip_tif(sb::SheetBuilder) = unzip_tif(full_folder_path(sb))
function unzip_tif(fofo)
    if isfile(joinpath(fofo, CONSOLIDATED_FNAM))
        @debug "$CONSOLIDATED_FNAM in $fofo already exists. Exiting `unzip_tif`."
        return true
    end
    zipfiles = filter(f-> endswith(f, ".zip"), readdir(fofo, join = true))
    for z in zipfiles
        _unzip_tif(z)
    end
    true
end

function _unzip_tif(zipfile)
    # Copied and modified for just .tif & .tiff from sylvaticus:
    # https://discourse.julialang.org/t/how-to-extract-a-file-in-a-zip-archive-without-using-os-specific-tools/34585/5
    @assert isfile(zipfile)
    basePath = dirname(zipfile)
    outPath = basePath
    @assert isdir(outPath)
    zarchive = ZipFile.Reader(zipfile)
    for f in zarchive.files
        fullFilePath = joinpath(outPath,f.name)
        if (endswith(f.name, ".tif") || endswith(f.name, ".tiff"))
            folder = dirname(fullFilePath)
            if ! ispath(folder)
                mkpath(folder)
            end
            if ! isfile(fullFilePath)
                open(fullFilePath, "w") do io
                    write(io, read(f))
                end
            end
        end
    end
    close(zarchive)
end


"""   
    tif_full_filenames_buried_in_folder(pth)

Excludes .tif files direcltly under 'metadata' folders
"""
function tif_full_filenames_buried_in_folder(pth)
    tifs = String[]
    for (root, dirs, files) in walkdir(pth)
        if ! endswith(root, "metadata")
            for file in files
                if endswith(file, ".tif")
                    fullname = joinpath(root, file)
                    @assert isfile(fullname)
                    push!(tifs, fullname)
                end
            end
        end
    end
    tifs
end