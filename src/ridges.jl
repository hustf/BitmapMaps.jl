# This is an experiment into a fluid dynamics analogy of terrain.

"""
    ridge_overlay(sb::SheetBuilder)
    ridge_overlay(fofo)

"""
function ridge_overlay(sb::SheetBuilder)
    ridge_overlay(full_folder_path(sb), sb.cell_iter, cell_to_utm_factor(sb), sb.density_pt_m⁻¹)
end
function ridge_overlay(fofo, cell_iter, cell2utm, density_pt_m⁻¹)
    if isfile(joinpath(fofo, RIDGE_FNAM))
        @debug "$RIDGE_FNAM in $fofo already exists. Exiting `ridge_overlay`."
        return true
     end
    if ! isfile(joinpath(fofo, CONSOLIDATED_FNAM))
        @debug "$CONSOLIDATED_FNAM in $fofo does not exist. Exiting `ridge_overlay`."
        return false
    end
    res = _ridge(fofo, cell_iter, cell2utm, density_pt_m⁻¹)
    ffna = joinpath(fofo, RIDGE_FNAM)
    @debug "Saving $ffna"
    save_png_with_phys(ffna, res, density_pt_m⁻¹)
    true
end
function _ridge(fofo, cell_iter, cell2utm, density_pt_m⁻¹)
    # Get elevation matrix. This samples every point regardless of cell_to_utm_factor
    g = readclose(joinpath(fofo, CONSOLIDATED_FNAM))
    za = g.A[:, :, 1]
    # Smooth surface
    za = imfilter(za, Kernel.gaussian(5));
    # Establish output image matrix
    ny, nx = size(cell_iter)
    source_indices = CartesianIndices((1:cell2utm:(nx * cell2utm), 1:cell2utm:(ny * cell2utm)))
    @debug "Render ridge"
    z_norm = za ./ cell2utm
    g1, g2 = imgradients(z_norm, KernelFactors.prewitt)
    g11, g12 = imgradients(g1, KernelFactors.prewitt)
    g21, g22 = imgradients(g2, KernelFactors.prewitt)
    divergence = g11 .+ g22
    dieder = RGBA{N0f8}(0.1, 0.0, 0.4, 0.5)
    corner = RGBA{N0f8}(0.0, 0.1, 0.1, 0.8)
    transpcol = RGBA{N0f8}(0, 0, 1, 0)
    ridge = map(source_indices) do I 
        x = divergence[I]
        if  x > 1 #0.1
            dieder
        elseif x < -1 # -0.05
            corner
        else
            transpcol
        end
    end
    # TODO: Remove small regions of dieders and ridges.
    #display(transpose(ridge))
    transpose(ridge)
end
