# Step in pipeline.
# Merges foreground layers onto background using the α / opacity channel.
# Output is an image file per sheet.

"""
    join_layers(sb::SheetBuilder)
    join_layers(fofo)

"""
join_layers(sb::SheetBuilder) = join_layers(full_folder_path(sb), sb.density_pt_m⁻¹)
function join_layers(fofo, density_pt_m⁻¹)
    ffna = joinpath(homedir(), fofo, COMPOSITE_FNAM)
    layerstack  = [TOPORELIEF_FNAM    Nothing
                   CONTOUR_FNAM      CompositeDestinationOver
                   RIDGE_FNAM        BlendMultiply           # This should be done in the HSL space to make ridges stand out more on snow.
                   WATER_FNAM        CompositeDestinationOver
                   GRID_FNAM         CompositeDestinationOver
                   MARKERS_FNAM      CompositeDestinationOver]
    layerfnas = joinpath.(homedir(), fofo, layerstack[:, 1])
    if isfile(ffna)
        # Do all layer files exist?
        if ! all(isfile.(layerfnas))
            @debug "    $COMPOSITE_FNAM in $fofo\n           already exists (not all layer files, though). Exiting `join_layers`"
            return true # We should continue with other sheets, no prob here
        end
        layermodtimes = mtime.(layerfnas)
        compmodtime = mtime(ffna)
        if ! any(layermodtimes .> compmodtime)
            @debug "    $COMPOSITE_FNAM in $fofo\n           already exists and is newer than layers. Exiting `join_layers`"
            return true
        end
    end
    if ! all(isfile.(layerfnas))
        @debug "    Layers missing. Exiting `join_layers`"
        return false
    end
    res = load(joinpath(homedir(), fofo, layerstack[1, 1]))
    # Iterate through the rest of the layer stack
    for rw in 2:size(layerstack, 1)
        layer = load(joinpath(homedir(), fofo, layerstack[rw, 1]))
        modefunc = layerstack[rw, 2]
        composite_on_top!(res, layer, modefunc)
    end
    # Feedback
    display_if_vscode(res)
    # Save file
    @debug "    Saving $ffna"
    save_png_with_phys(ffna, res; density_pt_m⁻¹)
    true
end
function composite_on_top!(res, layer, modefunc)
    @assert size(res) == size(layer)
    R = CartesianIndices(res)
    res[R] = modefunc.(layer[R], res[R])
    res
end

# TODO check out doing the layer multiplication in HSL colorspace.
# Currently, the ridge lines dominate forest.
#=
function composite_on_top!(res::AbstractArray{T1}, layer::AbstractArray{T1}, modefunc, colorspace::Type{T2}) where {T1, T2}
    @assert size(res) == size(layer)

    R = CartesianIndices(res)

    if T1 == T2
        # No color space conversion needed
        res[R] = modefunc.(layer[R], res[R])
    else
        # Perform blending with on-the-fly conversion
        res[R] = convert(T1, modefunc.(convert(T2, layer[R]), convert(T2, res[R])))
    end

    return res
end
=#