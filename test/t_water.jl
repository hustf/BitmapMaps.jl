# The algorithm and parameter values were optimized for a data set around 62.14111170231884N, 5.572971578864545E.
# Further refinements are hardcoded in this package.
#
# Work with the limited data in /resource. This area contains one tiny lake, hard to identify due
# to noise, but without some typical artefacts. A small corner of 'Holmevatnet' / Gurskøy is included,
# but the area is below the limit set for detection.

# Move the .tif file from /resource to a temporary directory.
using BitmapMaps
using Test

tmpdir_water = mktempdir()
let
    zfi = "../resource/eksport_800912_20240512.zip"
    zipfi = abspath(joinpath(@__DIR__, zfi))
    dest = joinpath(tmpdir_water, splitdir(zipfi)[2])
    cp(zipfi, dest)
end

# A zip file now exists in tmpdir_water, as if downloaded by user.
# Extract and inspect.
unzip_tif(tmpdir_water)
fna = first(tif_full_filenames_buried_in_folder(tmpdir_water))
# The file is zero padded with 1 m on the south boundary
@test polygon_string(fna) == "POLYGON ((18294 6937562, 18449 6937562, 18449 6937717, 18294 6937717, 18294 6937562),\n                   (18294 6937563, 18449 6937563, 18449 6937717, 18294 6937717, 18294 6937563))"

# This file is not 'consolidated', it is zero-valued on the southern edge. But the point here is
# not to test the pipeline, so we test on the level inside `water_overlay` first.
elevations = let
    z = readclose(fna)
    transpose(z.A[:, :, 1])
end
@test maximum(elevations) == 553.9368f0
@test eltype(elevations) == Float32

@test sum(iszero.(elevations)) == 155 # Southern edge zero

cell2utm = 1
cell_iter = CartesianIndices((1:cell2utm:155, 1:cell2utm:155))
lm_bool = BitmapMaps.is_water_surface(elevations, cell_iter, cell2utm, 0.155)
@test sum(lm_bool) * cell2utm^2 == 5372
cell2utm = 2
cell_iter = CartesianIndices((1:cell2utm:155, 1:cell2utm:155))
lm_bool = BitmapMaps.is_water_surface(elevations, cell_iter, cell2utm, 0.155)
@test sum(lm_bool) * cell2utm^2 == 5440 # Coarser along edges, OK

img = BitmapMaps.save_lakes_overlay_png(lm_bool, elevations, cell_iter, 1000, tmpdir_water)
@test isfile(joinpath(tmpdir_water, BitmapMaps.WATER_FNAM))
# The output file looks good. Now cleanup, then do the same through the 'pipeline_interface'
rm(joinpath(tmpdir_water, BitmapMaps.WATER_FNAM))

# Copy the elevation file to the filename expected by the pipeline interface.
cp(fna, joinpath(tmpdir_water, BitmapMaps.CONSOLIDATED_FNAM))
# A sheet builder fitting the data exactly
sb = let
    sheet_lower_left_utm = (18294, 6937562)
    pixel_origin_ref_to_bitmapmap = (0, 0)
    sheet_width_cell = 18449 -18294
    sheet_height_cell = 6937717 - 6937562
    cell_iter = CartesianIndices((1:sheet_height_cell, 1:sheet_width_cell))
    pthsh = tmpdir_water
    sheet_number = 1
    cell2utm = 1
    sheet_width_mm = 191
    density_pt_m⁻¹ = Int(ceil(sheet_width_cell  / sheet_width_mm))
    SheetBuilder(pixel_origin_ref_to_bitmapmap, cell_iter, sheet_lower_left_utm, cell2utm, sheet_number, density_pt_m⁻¹, pthsh)
end
# Pass this sheet builder to the interface, see if that makes us another WATER_FNAM file...
BitmapMaps.water_overlay(sb)
@test isfile(joinpath(tmpdir_water, BitmapMaps.WATER_FNAM))
