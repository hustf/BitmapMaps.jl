# Contains two steps in pipeline.
#
# 1. join_layers
# Merges foreground layers onto background using the α / opacity channel.
# Output is an image file per sheet.
#
# 2. make_thumbnail_image
# Unlike other pipeline functions, this takes a second, captured argument.

"""
    join_layers(sb::SheetBuilder)
    join_layers(fofo)

"""
join_layers(sb::SheetBuilder) = join_layers(full_folder_path(sb), sb.density_pt_m⁻¹)
function join_layers(fofo, density_pt_m⁻¹)
    if ! make_composite_image(fofo, density_pt_m⁻¹)
        return false
    end
    true
end
function make_composite_image(fofo, density_pt_m⁻¹)
    # Prepare 
    ffna = joinpath(homedir(), fofo, COMPOSITE_FNAM)
    layerstack  = [TOPORELIEF_FNAM    Nothing
                    CONTOUR_FNAM      CompositeDestinationOver
                    RIDGE_FNAM        BlendMultiply
                    WATER_FNAM        CompositeDestinationOver
                    GRID_FNAM         CompositeDestinationOver
                    MARKERS_FNAM      CompositeDestinationOver]
    layerfnas = joinpath.(homedir(), fofo, layerstack[:, 1])
    # Early exit
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
    for fna in layerfnas
        if ! isfile(fna)
            @debug "    Layers missing. Exiting `join_layers`: $(splitpath(fna)[end])"
            return false
        end
    end
    # Do the work
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

# This is called by by the generated function `make_thumbnail`, which captures n_governing
function make_thumbnail_image(fofo, density_pt_m⁻¹, n_governing)
    # Prepare 
    ffna = joinpath(homedir(), fofo, THUMBNAIL_FNAM)
    layerstack  = [TOPORELIEF_FNAM    Nothing
                    RIDGE_FNAM        BlendMultiply
                    WATER_FNAM        CompositeDestinationOver
                    MARKERS_FNAM      CompositeDestinationOver]
    layerfnas = joinpath.(homedir(), fofo, layerstack[:, 1])
    # Early exit
    if isfile(ffna)
        # Do all layer files exist?
        if ! all(isfile.(layerfnas))
            @debug "    $THUMBNAIL_FNAM in $fofo\n           already exists (not all layer files, though). Exiting `make_thumbnail`"
            return true # We should continue with other sheets, no prob here
        end
        layermodtimes = mtime.(layerfnas)
        compmodtime = mtime(ffna)
        if ! any(layermodtimes .> compmodtime)
            @debug "    $THUMBNAIL_FNAM in $fofo\n           already exists and is newer than layers. Exiting `make_thumbnail`"
            return true
        end
    end
    for fna in layerfnas
        if ! isfile(fna)
            @debug "    Layers missing. Exiting `make_thumbnail`: $(splitpath(fna)[end])"
            return false
        end
    end
    # Do the work
    img = load(joinpath(homedir(), fofo, layerstack[1, 1]))
    # Iterate through the rest of the layer stack
    for rw in 2:size(layerstack, 1)
        layer = load(joinpath(homedir(), fofo, layerstack[rw, 1]))
        modefunc = layerstack[rw, 2]
        composite_on_top!(img, layer, modefunc)
    end
    # Downsize
    ny, nx = size(img)
    iter = (1:n_governing:ny, 1:n_governing:nx)
    # A pure downsampling would result in a 'speckled' appearance. 
    # Averaging over the influencing pixels works well in the XYZ colorspace
    odd_size = 2 * div(n_governing, 2) + 1
    thumb = RGB{N0f8}.(mapwindow(mean, XYZ.(img), (odd_size, odd_size); indices = iter))
    # Feedback
    display_if_vscode(thumb)
    # Save file
    @debug "    Saving $ffna"
    save_png_with_phys(ffna, thumb; density_pt_m⁻¹)
    true
end
