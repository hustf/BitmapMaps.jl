# The algorithm and parameter values is optimized for in environments/tutorial_images/image_segmentation.jl
# Since these parameter values are well optimized, we use those hardcoded here.
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
    zipfi = joinpath(@__DIR__, zfi)
    cp(zipfi, joinpath(tmpdir_water, splitdir(zipfi)[2]))
end

# A zip file now exists in tmpdir_water, as if downloaded by user.
# Extract and inspect.
unzip_tif(tmpdir_water) 
fna = first(tif_full_filenames_buried_in_folder(tmpdir_water))
@test BitmapMaps.bounding_box_polygon_string(fna)== "POLYGON ((18294 6937562, 18449 6937562, 18449 6937717, 18294 6937717, 18294 6937562))"

# This file is already 'consolidated', i.e. dense. The point here is not to test the pipeline, so we 
# test on the level inside `water_overlay` first.

elevations = let 
    z = BitmapMaps.GeoArrays.read(fna)
    transpose(z.A[:, :, 1])
end
@test maximum(elevations) == 553.9368f0
@test eltype(elevations) == Float32
lm_bool = BitmapMaps.is_water_surface(elevations, 1)
@test sum(lm_bool) == 5205 # This lake's identified surface are is 5205m²
img = BitmapMaps.save_lakes_overlay_png(lm_bool, elevations, 1000, tmpdir_water)
@test isfile(joinpath(tmpdir_water, BitmapMaps.WATER_FNAM))
# The output file looks good. Now cleanup, then do the same through the 'pipeline_interface'
rm(joinpath(tmpdir_water, BitmapMaps.WATER_FNAM))

# Copy the elevation file to the filename expected by the pipeline interface.
# (we would move the file, but unfortunately it is not closed by the C library we use) 
cp(fna, joinpath(tmpdir_water, BitmapMaps.CONSOLIDATED_FNAM))
# A sheet builder fitting the data exactly
sb = let
    sheet_lower_left_utm = (18294, 6937562)
    pixel_origin_ref_to_bitmapmap = (0, 0)
    sheet_pix_width = 18449 -18294
    sheet_pix_height = 6937717 - 6937562
    pix_iter = CartesianIndices((1:sheet_pix_height, 1:sheet_pix_width))
    pthsh = tmpdir_water
    sheet_number = 1
    pix_dist = 1
    SheetBuilder(sheet_lower_left_utm, pix_dist, pixel_origin_ref_to_bitmapmap, pix_iter, pthsh, sheet_number)
end
# Pass this sheet builder to the interface, see if that makes us another WATER_FNAM file...
BitmapMaps.water_overlay(sb)
@test isfile(joinpath(tmpdir_water, BitmapMaps.WATER_FNAM))
