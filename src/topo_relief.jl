# This renders a topographic relief map.

#=
sheet_width_cell = floor(sheet_width_mm * density_pt_m⁻¹ / 25.4)
sheet_width_cell = floor(sheet_width_mm * density_pt_m⁻¹ / 25.4)
density_pt_m⁻¹ = density_pt_m⁻¹ / 0.0254
0.0254 * density_pt_m⁻¹ = density_pt_m⁻¹
sheet_width_cell = floor(sheet_width_mm * density_pt_m⁻¹ / 1000)
sheet_width_cell = floor(sheet_width_mm * density_pt_m⁻¹ / 1000)
=#

"""
    topo_relief(sb::SheetBuilder)
    topo_relief(fofo)

"""
topo_relief(sb::SheetBuilder) = topo_relief(full_folder_path(sb), sb.cell_iter, sb.sheet_width_mm)
function topo_relief(fofo, cell_iter, sheet_width_mm)
    if isfile(joinpath(fofo, TOPORELIEF_FNAM))
        @debug "$TOPORELIEF_FNAM in $fofo already exists. Exiting `topo_relief`."
        return true
     end
    if ! isfile(joinpath(fofo, CONSOLIDATED_FNAM))
        @debug "$CONSOLIDATED_FNAM in $fofo does not exist. Exiting `topo_relief`."
        return false
    end
    # TODO take the value as an argument instead.
    density_pt_m⁻¹ = Int(round(1000 * size(cell_iter)[2] / sheet_width_mm))
    _topo_relief(fofo, cell_iter, density_pt_m⁻¹)
end
function _topo_relief(fofo, cell_iter, density_pt_m⁻¹)
    # Get elevation matrix
    za = let 
        z = readclose(joinpath(fofo, CONSOLIDATED_FNAM))
        transpose(z.A[:, :, 1])
    end
    @assert maximum(za) > 0
    # Establish output image matrix
    rgbimg = zeros(RGBA{N0f8}, size(cell_iter))
    @assert size(za) == size(rgbimg)
    mi, ma = extrema(za)
    mi = 541
    #@show mi ma

    for I in cell_iter
        @assert za[I] >= 0
 #       zrel = N0f8(za[I] / 1500)
        zr = max(0.0, (za[I] - mi) / (ma - mi))
        zrel = N0f8(zr)
        if zr == 0
            rgbimg[I] = RGBA{N0f8}(1.0, 0.5, 0.0, 1.0)
        else
            rgbimg[I] = RGBA{N0f8}(zrel, zrel, zrel, 1.0)
        end
    end 
    save_png_with_phys(joinpath(fofo, BitmapMaps.TOPORELIEF_FNAM), rgbimg, density_pt_m⁻¹)
    true
end

function possible_integer_densities(sheet_width_cell, sheet_width_mm)
    throw("dead")
    # Initialize an array to hold possible values of density_pt_m⁻¹
    densities = Int[]

    # Calculate the minimum and maximum integer values for density_pt_m⁻¹ that could produce the given sheet_width_cell
    min_density = ceil(sheet_width_cell * 1000 / sheet_width_mm)
    max_density = floor((sheet_width_cell + 1) * 1000 / sheet_width_mm) - 1

    # Loop through the possible range and collect valid densities
    for density in min_density:max_density
        if floor(sheet_width_mm * density / 1000) == sheet_width_cell
            push!(densities, density)
        end
    end
    
    return densities
end