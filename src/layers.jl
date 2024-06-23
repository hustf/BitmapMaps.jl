# Merges foreground layers onto background using the α / opacity channel.
# Output is an image file per sheet.

"""
    join_layers(sb::SheetBuilder)
    join_layers(fofo)

"""
join_layers(sb::SheetBuilder) = join_layers(full_folder_path(sb), sb.cell_iter, sb.density_pt_m⁻¹)
function join_layers(fofo, cell_iter, density_pt_m⁻¹)
    if isfile(joinpath(fofo, COMPOSITE_FNAM))
        @debug "$COMPOSITE_FNAM in $fofo already exists. Exiting `join_layers`."
        return true
    end
    if ! isfile(joinpath(fofo, TOPORELIEF_FNAM))
        @debug "$TOPORELIEF_FNAM in $fofo does not exist. Exiting `join_layers`."
        return false
    end
    res = _join_layers(fofo, cell_iter)
    ffna = joinpath(fofo, COMPOSITE_FNAM)
    @debug "Saving $ffna"
    save_png_with_phys(ffna, res, density_pt_m⁻¹)
    true
end
function _join_layers(fofo, cell_iter)
    # Note, an in-place version would be faster.
    comp = composite_image(joinpath(fofo, TOPORELIEF_FNAM), joinpath(fofo, WATER_FNAM))
    comp = composite_image(comp, joinpath(fofo, CONTOUR_FNAM))
    comp = composite_image(comp, joinpath(fofo, GRID_FNAM))
    display(comp)
    return comp
end


function composite_image(fna_s::String, fna_d::String)
    s = load(fna_s)
    d = load(fna_d)
    composite_image(s, d)
end

function composite_image(s, fna_d::String)
    d = load(fna_d)
    composite_image(s, d)
end
function composite_image(s, d)
    if size(s) == size(d)
        map(CompositeDestinationOver, d, s)
    else
        ny, nx = size(s)
        # Will error if s is larger than d
        cell2utm = Int(size(d)[1] / ny)
        source_indices = CartesianIndices((1:cell2utm:(ny * cell2utm), 1:cell2utm:(nx * cell2utm)))
        shrunk_d = map(I -> d[I], source_indices)
        map(CompositeDestinationOver, shrunk_d, s)
    end
end
