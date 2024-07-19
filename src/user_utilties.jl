# The pipeline has a consolidation step, but the pipeline requires that user download relevant elevation data from online.
#
# However, if data already exists locally for another bitmapmap, these functions
# can help with local file management.
#
# To build the folder structure first, you can either `run_bitmapmap_pipeline()`,
# which aborts the run when not finding data, or run
# `establish_folder.(sheet_matrix_builder)`
"""
    copy_relevant_tifs_to_folder((source_folder, destination_folder)
    copy_relevant_tifs_to_folder(source_folder, smb::SheetMatrixBuilder)
    copy_relevant_tifs_to_folder(source_folder, sb::SheetBuilder)
    ---> Vector{String}, full names of new files in destination

- `source_folder` is searched recursively. Each file is opened to find if it's relevant to the destination coordinates.
- `destination_folder` respects the integers naming scheme: `r c min_x min_y max_x max_y` (ref. `parse_folder_name`).
   Relevant files copied in the first level of destination's hierarchy.

Any files to be copied need to match the geographical area specified in destination folder name.

Ignores any files named `Bitmapmap.CONSOLIDATED_FNAM`.

Copying from somewhere inside the destination folder to the destination folder is not possible, as this would duplicate files for no good reason.
A more advanced version would create shortcuts, but we don't.
"""
function copy_relevant_tifs_to_folder(source_folder, destination_folder)
    # Any files to be copied need to match the geographical area of source.
    r, c, min_x, min_y, max_x, max_y = parse_folder_name(destination_folder)
    d = (;min_x = Float64(min_x), min_y = Float64(min_y), max_x = Float64(max_x), max_y = Float64(max_y))
    # Candidates from file names only, no geographical info. If the same file name occurs twice in the source folder hierarchy,
    # only one will be copied.
    cs = candidate_tif_names(source_folder, destination_folder)
    @info "Found $(length(cs)) candidates for copying. Set ENV[\"JULIA_DEBUG\"] = \"BitmapMaps\" for detailed output."
    # Candidate by candidate is checked for geographical match, then copied.
    destination_files = String[]
    for cafna in cs
        printstyled(stdout, repeat(' ', 8), cafna , "\n"; color = :light_black)
        if is_source_relevant(d, cafna)
            @debug "    Found relevant file $cafna"
            dfna = file_name_at_dest(cafna, destination_folder)
            cp(cafna, dfna)
            push!(destination_files, dfna)
        else
            bbs = closed_box_string(d)
            @debug "    Ignored copying $(splitdir(cafna)[end]) because it does not overlap the geographical region $bbs"
        end
    end
    destination_files
end
function copy_relevant_tifs_to_folder(source_folder, smb::SheetMatrixBuilder)
    destination_files = String[]
    for sb in smb
        append!(destination_files, copy_relevant_tifs_to_folder(source_folder, sb))
    end
    destination_files
end
function copy_relevant_tifs_to_folder(source_folder, sb::SheetBuilder)
    copy_relevant_tifs_to_folder(source_folder, full_folder_path(sb))
end


function file_name_at_dest(full_file_name, destination_folder)
    fna_short = splitdir(full_file_name)[2]
    joinpath(destination_folder, fna_short)
end

"""
    candidate_tif_names(source_folder, destination_folder))

Candidates for copying based on file names, excepting file names that already exists on
any level of the destination folder hierarchy.
"""
function candidate_tif_names(source_folder, destination_folder)
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
    fnas_candidates = tif_full_filenames_buried_in_folder(source_fo)
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

function display_if_vscode(M::Matrix{T}) where T <: Union{RGBA{N0f8}, RGB{N0f8}}
    if isinteractive()
        if get(ENV, "TERM_PROGRAM", "") == "vscode"
            # Display
            display(M)
        end
    end
end