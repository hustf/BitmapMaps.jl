# The pipeline has a consolidation step, but the pipeline requires that user download relevant elevation data from online.
#
# However, if data already exists locally for another bitmapmap, these functions
# can help with local file management.
#
# To build the folder structure first, you can either `run_bitmapmap_pipeline()`,
# which aborts the run when not finding data, or run
# `establish_folder.(sheet_matrix_builder)`

"""
    copy_relevant_tifs_to_folder((source_folder, destination_folder; recurse = true)
    copy_relevant_tifs_to_folder(source_folder, smb::SheetMatrixBuilder; recurse = true)
    copy_relevant_tifs_to_folder(source_folder, sb::SheetBuilder; recurse = true)
    ---> Vector{String}, full names of new files in destination

- `source_folder` is searched recursively if 'recurse' is true. Each file is opened to find if it's relevant to the destination coordinates.
- `destination_folder` respects the integers naming scheme: `r c min_x min_y max_x max_y` (ref. `parse_folder_name`).
   Relevant files copied into the destination's first level.

Any files to be copied need to match the geographical area specified in destination folder name.

Ignores any files named `Bitmapmap.CONSOLIDATED_FNAM`.

Copying from somewhere inside the destination folder to the destination folder is not possible, as this would duplicate files for no good reason.
A more advanced version would create shortcuts, but we don't.
"""
function copy_relevant_tifs_to_folder(source_folder, destination_folder; recurse = true)
    # Any files to be copied need to match the geographical area of source.
    r, c, min_x, min_y, max_x, max_y = parse_folder_name(destination_folder)
    d = (;min_x = Float64(min_x), min_y = Float64(min_y), max_x = Float64(max_x), max_y = Float64(max_y))
    # Candidates from file names only, no geographical info. If the same file name occurs twice in the source folder hierarchy,
    # only one will be copied.
    cs = candidate_tif_names(source_folder, destination_folder; recurse)
    if get(ENV, "JULIA_DEBUG", "") !== "BitmapMaps"
        @info "Found $(length(cs)) candidates for copying."
    else
        @info "Found $(length(cs)) candidates for copying. Set ENV[\"JULIA_DEBUG\"] = \"BitmapMaps\" for detailed output."
    end
    # Candidate by candidate is checked for geographical match, then copied.
    destination_files = String[]
    for cafna in cs
        printstyled(stdout, repeat(' ', 8), cafna , "\n"; color = :light_black)
        if is_source_relevant(d, cafna)
            @debug "    Found relevant file $cafna"
            dfna = file_name_at_dest(cafna, destination_folder)
            if isfile(dfna) 
                @warn "Existing file: $dfna, hinders copying $cafna "
            else
                @debug "    Destination file $dfna"
                # This has at least once failed for a file > 1 Gb, but not every time. Do it manually if that occurs.
                cp(cafna, dfna)
                push!(destination_files, dfna)
            end
        else
            bbs = closed_box_string(d)
            @debug "    Ignored copying $(splitdir(cafna)[end]) because it does not overlap the geographical region $bbs"
        end
    end
    destination_files
end
function copy_relevant_tifs_to_folder(source_folder, smb::SheetMatrixBuilder; recurse = true)
    destination_files = String[]
    for sb in smb
        append!(destination_files, copy_relevant_tifs_to_folder(source_folder, sb; recurse))
    end
    destination_files
end
function copy_relevant_tifs_to_folder(source_folder, sb::SheetBuilder; recurse = true)
    copy_relevant_tifs_to_folder(source_folder, full_folder_path(sb); recurse)
end


function file_name_at_dest(full_file_name, destination_folder)
    fna_short = splitdir(full_file_name)[2]
    joinpath(destination_folder, fna_short)
end

"""
    candidate_tif_names(source_folder, destination_folder; recurse = true)

Candidates for copying based on file names, excepting file names that already exists on
any level of the destination folder hierarchy.

If 'recurse' is true, will look in the hierarchy beneath source_folder .
"""
function candidate_tif_names(source_folder, destination_folder; recurse = true)
    source_fo = abspath(source_folder)
    dest_fo = abspath(destination_folder)
    # For avoiding unwanted duplicates:
    fnas_avoid = let
        # These exist, but maybe buried in the destination folders.
        fnas_dest = tif_full_filenames_buried_in_folder(dest_fo)
        # We'd like to avoid making duplicates of those at the top level of destination_folder:
        map(fna-> file_name_at_dest(fna, dest_fo), fnas_dest)
    end
    # Also avoid the name of the file we're going to make later:
    push!(fnas_avoid, file_name_at_dest(BitmapMaps.CONSOLIDATED_FNAM, dest_fo))
    # These exist, but may be buried in the destination folders.
    fnas_candidates = tif_full_filenames_buried_in_folder(source_fo; recurse)
    filter(fnas_candidates) do cand
        file_name_at_dest(cand, dest_fo) ∉ fnas_avoid
    end
end



"""
    bbox_external_overlap(
        bbox_a::NamedTuple{(:min_x, :min_y, :max_x, :max_y)},
        bbox_b::NamedTuple{(:min_x, :min_y, :max_x, :max_y)})
        ---> Bool

This assumes the boundary boxes are 'external', i.e. enclosing cells that are part of it.

Returns false for adjacent bboxes, like e.g.
```
julia> BitmapMaps.bbox_external_overlap(
    (min_x = 44001, min_y = 47, max_x = 44006, max_y = 55),
    (min_x = 44006, min_y = 47, max_x = 44012, max_y = 55))
false
````
"""
function bbox_external_overlap(
    a::NamedTuple{(:min_x, :min_y, :max_x, :max_y)},
    b::NamedTuple{(:min_x, :min_y, :max_x, :max_y)})
    # Check if external bboxes are valid
    if (a.min_x > a.max_x) || (a.min_y > a.max_y)
        throw(ArgumentError("Invalid a (min > max) $a"))
    end
    if (b.min_x > b.max_x) || (b.min_y > b.max_y)
        throw(ArgumentError(("Invalid b (min > max) $b")))
    end
    # Check if bboxes do not overlap
    if (a.max_x <= b.min_x) ||
       (a.max_y <= b.min_y) ||
       (a.min_x >= b.max_x) ||
       (a.min_y >= b.max_y)
        return false
    end
    true
end

function display_if_vscode(M)
    if isinteractive()
        if get(ENV, "TERM_PROGRAM", "") == "vscode"
            # Stretch gray colors from black to white
            foo = scaleminmax(extrema(M)...)
            # Display
            display(colorview(Gray, foo.(M)))
        end
    end
end

function display_if_vscode(M::Matrix{T}) where T <: Union{RGBA{N0f8}, RGB{N0f8}, XYZ{Float32}, XYZA{Float32}}
    if isinteractive()
        if get(ENV, "TERM_PROGRAM", "") == "vscode"
            # Display
            display(M)
        end
    end
end



"""
    write_vectors_to_csv(filename::String, headers::Vector{T}, vectors, widths; delim = '\t') where T <: Union{String, Symbol}
    ---> Nothing

Make a column formatted and delimited .csv, without many dependencies. Not suitable for large files.

# Example
```
julia> write_vectors_to_csv(filename::String, headers::Vector{String}, vectors, widths)
```
"""
function write_vectors_to_csv(filename::String, headers::Vector{T}, vectors, widths; delim = '\t') where T <: Union{String, Symbol}
    length(headers) == length(vectors) == length(widths) ||  throw(ArgumentError("NOT length(vectors) == length(headers) == length(widths)"))
    # Check that all vectors are of the same length
    length_of_vectors = length(vectors[1])
    for v in vectors
        if length(v) != length_of_vectors
            throw(ArgumentError("All vectors must be of the same length"))
        end
    end
    sh = join(rpad.(string.(headers), widths), delim)
    open(filename, "w") do io
        # Write header
        write(io, sh * "\n")
        # Write data
        for i in 1:length_of_vectors
            row = [vectors[j][i] for j in 1:length(vectors)]
            s = join(rpad.(string.(row), widths), delim)
            write(io, s * "\n")
        end
    end
    nothing
end


# Function useful for developing and testing

function open_as_temp__in_gimp(img)
    @async let
        fnam = tempname()
        save_png_with_phys(fnam, img)
        gimp_path = "C:\\Program Files\\GIMP 2\\bin\\gimp-2.10.exe"
        run(`$gimp_path "$fnam"`)
    end
end

function get_random_color(i)
    Random.seed!(i) # For consistentency between runs
    rand(RGB{N0f8})
end