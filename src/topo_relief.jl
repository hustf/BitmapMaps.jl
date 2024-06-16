# This renders a topographic relief map.

"""
    topo_relief(sb::SheetBuilder)
    topo_relief(fofo)

"""
topo_relief(sb::SheetBuilder) = topo_relief(full_folder_path(sb), sb.cell_iter, sb.density_pt_m⁻¹)
function topo_relief(fofo, cell_iter, density_pt_m⁻¹)
    if isfile(joinpath(fofo, TOPORELIEF_FNAM))
        @debug "$TOPORELIEF_FNAM in $fofo already exists. Exiting `topo_relief`."
        return true
     end
    if ! isfile(joinpath(fofo, CONSOLIDATED_FNAM))
        @debug "$CONSOLIDATED_FNAM in $fofo does not exist. Exiting `topo_relief`."
        return false
    end
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

