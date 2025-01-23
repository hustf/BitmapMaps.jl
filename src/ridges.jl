# Step in pipeline.
# Calculates the divergence in a fluid flow analogy, where the
# gradient of elevation resembles fluid flow. The same operator 
# is called 'laplacian'.
#
# Large divergence coincicides with terrain ridges or outward corners.
# Large negative divergence coincides with ridges, canyons or inward corners.
#
# In built-up and foresty areas, this is masked by local roughness, 
# so we apply the same mask as is used in 'countours':  bumpy_patch.
# The ridge and dieder lines are not generated in bumpy patches.
#
# Output is an image file per sheet, for manual touch-up.
# Consider TODO: Many artifacts have a very even, small, positive curvature.
# We might identify those in a similar way as we do with lakes? 

"""
    ridge_overlay(sb::SheetBuilder)
    ridge_overlay(fofo)

"""
function ridge_overlay(sb::SheetBuilder)
    ridge_overlay(full_folder_path(sb), sb.cell_iter, cell_to_utm_factor(sb))
end
function ridge_overlay(fofo, cell_iter, cell2utm)
    if isfile(joinpath(fofo, RIDGE_FNAM))
        @debug "    $RIDGE_FNAM in $fofo \n           already exists. Exiting `ridge_overlay`"
        return true
     end
    if ! isfile(joinpath(fofo, CONSOLIDATED_FNAM))
        @debug "    $CONSOLIDATED_FNAM in $fofo\n           does not exist. Exiting `ridge_overlay`"
        return false
    end
    res = _ridge(fofo, cell_iter, cell2utm)
    # Feedback
    display_if_vscode(res)
    # Save
    ffna = joinpath(fofo, RIDGE_FNAM)
    @debug "    Saving $ffna"
    save_png_with_phys(ffna, res)
    true
end
function _ridge(fofo, cell_iter, cell2utm)
    # Get elevation matrix. This samples every point regardless of cell_to_utm_factor
    g = readclose(joinpath(fofo, CONSOLIDATED_FNAM))
    @debug "    Render ridge"
    # Define output size and colours
    ny, nx = size(cell_iter)
    source_indices = (1:cell2utm:(ny  * cell2utm), 1:cell2utm:(nx * cell2utm))
    si = CartesianIndices(source_indices)
    # Find the scalar matrix representing 'negative divergence'.
    # Interpretation: Gradients points uphill.
    # If the terrain gradients were a 2d fluid vector field, 
    # a negative laplacian at a point means that
    # 'fluid' collects at a point (the divergence is positive).
    #
    # Where there is forest or houses, we do not want to show ridges.
    # masked_laplacian will be zero in those patches.
    # TODO consider using σₘᵢₙ , σₘₐₓ = principal_min_max_stress(z), see 'dev_summit_markers.jl'
    masked_laplacian = smooth_laplacian(g, si) .* (1 .- bumpy_patch(g, si))
    # Pre-allocate output image
    result = zeros(RGBA{N0f8}, size(si)...)
    # Pre-allocate boolean buffer. This has the same size as the output image
    bbuf = Array{Gray{Bool}}(undef, size(si)...)
    # Add lines to 'result', different colors for different criteria.
    let 
        colo_corner = RGBA{N0f8}(0.118, 0.102, 0.141, 0.85) # WAS RGBA{N0f8}(tuple(easter_map_shadow_colors()[end])..., 0.85)
        colo_dieder = RGBA{N0f8}(0.0, 0.64, 1.0, 0.7)  # WAS (0.02, 0.0, 0.075, 0.7) # WAS (tuple(easter_map_shadow_colors()[1])..., 0.7)
        criterion_functions = [<(-0.05f0), >(0.2f0)]
        thicknesses = [3, 5]
        _corners_and_dieders!(result, bbuf, masked_laplacian, criterion_functions, [colo_corner, colo_dieder], thicknesses)
    end
    result
end
function _corners_and_dieders!(result, bbuf, source, criterion_functions, colors, thicknesses)
    @assert length(criterion_functions) == length(colors) == length(thicknesses) == 2
    for (colo, t, criterion_function) in zip(colors, thicknesses, criterion_functions)
        # Overwrite bbuf with pixels on corners
        mapwindow!(M -> Gray{Bool}(criterion_function(first(M))), bbuf, source, (1, 1))
        # Our very clever thinning, despeckling and thickening:
        bbuf .= strokeify(bbuf, t, 5) # Minimum length of a line is hardcoded as 5.
        # Edges of sheets sometimes get stroked due to boundary conditions. 
        # It's better to clean the edges fully.
        bbuf[:, 1:t] .= false
        bbuf[1:t, :] .= false
        bbuf[(end - t):end, :] .= false
        bbuf[:, (end - t):end] .= false
        # Treat 'false' as transparent. Overlay bbuf on result and write in-place to result.
        map!((outcol, bol) -> bol == true ? colo : outcol, result, result, bbuf)
    end
    result
end

"""
    σm(z::Matrix{Float32})
    ---> Matrix{Float32}

Mean stress, i.e. hydrostatic pressure with the opposite sign. Tensile stress
is positive.

Use this for identifying cone-like summits, conifers or towers which have negative values.
Use this for identifying holes in the ground or bowls, which have positive values.
"""
function σm(z::Matrix{Float32})
    # The Hessian components are used to determine the principal tensile stress equivalent
    σ11, _, _, σ22 = hessian_components(z)
    # Return a scalar matrix representing the mean stress at each index
    (σ11 .+ σ22) ./ 2
end


function principal_stresses(σ11, σ12, σ22)
    trace = (σ11 .+ σ22) ./ 2
    determinant = sqrt.(((σ11 .- σ22) ./ 2).^2 .+ σ12.^2)
    λ1 = trace .+ determinant
    λ2 = trace .- determinant
    λ1, λ2
end
function principal_stresses(z)
    σ11, _, σ12, σ22 = hessian_components(z)
    # Disregarding roundoff, σ21 == σ22
    principal_stresses(σ11, σ12, σ22)
end

principal_stress_min(z) = min.(principal_stresses(z))
principal_stress_max(z) = max.(principal_stresses(z))

function principal_min_max_stress(z)
    λ1, λ2 = principal_stresses(z)
    min.(λ1, λ2), max.(λ1, λ2)
end




"""
   smooth_laplacian(g, source_indices::CartesianIndices)
   ---> Matrix{scalar}

We can interpret the gradient of elevation z as a vector field.

If that vector field represented a quasi-static, incompressible
fluid flow, we would call the result 'divergence': How much
fluid would be entering (from outside) at this point? 
"""
function smooth_laplacian(g, source_indices::CartesianIndices)
    # Smooth surface
    zf = smooth_surface_fir(g; w = 159, nyquist_denom = 16)
    g11, _, _, g22 = hessian_components(zf)

    # We have used the full resolution to determine the curvature.
    # Return the needed resolution only.
    g11[source_indices] .+ g22[source_indices]
end

jacobian_components(z) = imgradients(z, KernelFactors.prewitt)

"""
    hessian_components(g1, g2)
    ---> g11, g21, g12, g22 {Float32}

g1, g2 are intended to be matrices containing orthogonal gradient components.

Output is matrices of the same form, representing second order derivatives.
These are based on the Prewitt algorithm, i.e. they are 'smoothed sideways'.

Each point [0,0] is affected by region [-2:2, 2:2].
"""
function hessian_components(g1, g2)
    # Hessian matrix elements (g12 and g21 are supposedly equal)
    g11, g12 = imgradients(g1, KernelFactors.prewitt)
    g21, g22 = imgradients(g2, KernelFactors.prewitt)
    g11, g21, g12, g22
end
function hessian_components(z::Matrix{T}) where T
    # Jacobian matrix elements
    g1, g2 = jacobian_components((z))
    hessian_components(g1, g2)
end
