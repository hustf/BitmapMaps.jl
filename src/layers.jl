# Contains two steps in pipeline. Also: `downsize` and `downsample`.
#
# 1. Join_layers
# Merges foreground layers onto background using the α / opacity channel.
# Output is an image file per sheet.
#
# 2. Make_thumbnail_image
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

"""
    make_thumbnail_image(fofo, density_pt_m⁻¹, n_governing)
    ---> Bool

Called by by the generated function `make_thumbnail`, which captures n_governing.

The thumbnail contains fewer detail layers, and uses downsizing instead of downsampling.
"""
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
        # We will not update the thumbnail just because MARKERS_FNAM is newer than the composite image.
        # The reason why is that 'MARKERS_FNAM' will be updated in 'summits_regional_poststep'.
        layermodtimes = mtime.(layerfnas[1:end - 1])
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
    thumb = downsize(img, n_governing)
    # Feedback
    display_if_vscode(thumb)
    # Save file
    @debug "    Saving $ffna"
    save_png_with_phys(ffna, thumb; density_pt_m⁻¹)
    true
end

"""
    downsize(img::Matrix{T}, n_governing) where T
    downsize(img::Matrix{T}, n_governing) where T<:Union{Float32, Float64}
    ---> typeof(img)

Downsizing returns the mean of the original pixels for an output pixel.

For colors, it averages in the XYZ colorspace, since linear operations in RGB do not
match visual perception.

"""
function downsize(img::Matrix{T}, n_governing) where T
    ny, nx = size(img)
    iter = (1:n_governing:ny, 1:n_governing:nx)
    odd_size = 2 * div(n_governing, 2) + 1
    T.(mapwindow(mean, XYZ.(img), (odd_size, odd_size); indices = iter))
end
function downsize(img::Matrix{T}, n_governing) where T<:Union{Float32, Float64}
    ny, nx = size(img)
    iter = (1:n_governing:ny, 1:n_governing:nx)
    # A pure downsampling would result in a 'speckled' appearance.
    # Averaging over the influencing pixels
    odd_size = 2 * div(n_governing, 2) + 1
    mapwindow(mean, img, (odd_size, odd_size); indices = iter)
end

"""
    downsample(img, n_governing)

Downsampling retains sharpness and is preferrable to
`downsize` when the result is going to be compared to a threshold.
"""
function downsample(img, n_governing)
    ny, nx = size(img)
    iter = CartesianIndices((1:n_governing:ny, 1:n_governing:nx))
    img[iter]
end