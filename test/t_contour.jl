# Internals to _elev_contour.
# This is best checked in a slightly large area, which we won't 
# include in the repo resources. Instead, we rely on local files,
# which can be downloaded without creating a user account.
using Test
using BitmapMaps
using BitmapMaps: fir_lp_coefficients, centered
using BitmapMaps: imfilter, mapwindow, mapwindow!, define_builder, _elev_contours
using BitmapMaps: contour_lines_overlay, CONTOUR_FNAM

#fofo = joinpath(homedir(), "bitmapmaps\\render\\1 1  47675 6929520  54041 6938686")
#fofo = joinpath(homedir(), raw"BitmapMaps/Hov Litlevatn\1 1  29000 6826535  31828 6830608")
fofo = joinpath(homedir(), raw"BitmapMaps\proj 47675 6929520 57224 6947852\1 1  35425 6920995  42190 6930739")
if ! ispath(joinpath(fofo))
    throw("Missing folder, cannot continue. Try to establish $fofo and \n copy_relevant_tifs_to_folder(source_folder, destination_folder) ")
end
if ! isfile(joinpath(fofo, BitmapMaps.CONSOLIDATED_FNAM))
    throw("Missing file, cannot continue. Try \n BitmapMaps.consolidate_local_data_to_geoarray_in_folder(fofo)")
end
minlen = 6
vthick = [1, 3, 5]
vdist = [20, 100, 1000]
cell2utm = 3
nx = 2255
ny = 3248
cell_iter = CartesianIndices((1:ny, 1:nx))

# With v0.1.3: 11.523156 seconds (67.75 k allocations: 9.224 GiB, 16.94% gc time)
# With func_I_to_utm:  11.518735 seconds (67.75 k allocations: 9.224 GiB, 16.91% gc time)
# After reducing resolution, temporary state: 8.394829 seconds (15.79 k allocations: 5.546 GiB, 6.76% gc time)
# Current: 9.757608 seconds (22.45 M allocations: 7.385 GiB, 8.24% gc time)
@time res = _elev_contours(fofo, cell_iter, cell2utm, minlen, vthick, vdist)

rm(joinpath(joinpath(fofo, CONTOUR_FNAM)))
@time contour_lines_overlay(fofo, cell_iter, cell2utm, minlen, vthick, vdist)















