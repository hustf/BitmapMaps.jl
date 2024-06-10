using Test
using BitmapMaps

##################
# user_utilties.jl
##################

tempdir = mktempdir()
cd(tempdir)
source_folder = "s"
mkpath(source_folder)
mkpath("$source_folder/a")
mkpath("$source_folder/b")
open("$source_folder/a/1.tif", "w") do io
    println(io, "A tif text file")
end
open("$source_folder/b/2.tif", "w") do io
    println(io, "A tif text file")
end
open("$source_folder/b/$(BitmapMaps.CONSOLIDATED_FNAM)", "w") do io # This file name is reserved, won't be copied.
    println(io, "A tif text file")
end
# Dest 
dfona = "2 3 100 1000 200 2000"
mkpath(dfona)
@test length(BitmapMaps.candidate_tif_names(source_folder, dfona)) == 2

# These "text file tif"s didn't have the right format, though:
@test_throws Exception copy_relevant_tifs_to_folder(joinpath(tempdir, source_folder), dfona)

# Remove the fake tifs
rm("$source_folder/a/1.tif")
rm("$source_folder/b/2.tif")
# Let's unpack a zip files from the /resource folder, which should be relevant.
zipfi = joinpath(@__DIR__, "../resource/eksport_796345_20240420.zip")
pth = "$source_folder/a/"
fzip = joinpath(pth, splitdir(zipfi)[2])
cp(zipfi, fzip)
# A downloaded and relevant .zip file now exists in source_folder/a. Unzip it to two  .tif files
unzip_tif(pth)
@test length(BitmapMaps.candidate_tif_names(source_folder, dfona)) == 2
# But these .tif files cover another area than defined by our folder name
copy_relevant_tifs_to_folder(joinpath(tempdir, source_folder), dfona)
@test length(BitmapMaps.candidate_tif_names(source_folder, dfona)) == 2

# This folder name overlaps one of the tif files.
dfona = "5 6 43300 6909000 44000 6909200"
mkpath(dfona)
@test length(BitmapMaps.candidate_tif_names(source_folder, dfona)) == 2
copy_relevant_tifs_to_folder(joinpath(tempdir, source_folder), dfona)
@test length(BitmapMaps.candidate_tif_names(source_folder, dfona)) == 2

# Methods     
#    copy_relevant_tifs_to_folder(source_folder, smb::SheetMatrixBuilder)
#    copy_relevant_tifs_to_folder(source_folder, sb::SheetBuilder)
# are tested in 't_pipeline.jl'

######################
# geoarray_utilties.jl
######################

fnas = tif_full_filenames_buried_in_folder(pth)
@test nonzero_raster_closed_polygon_string(fnas[1]) == "(44000 6909047, 44001 6909047, 44001 6909055, 44000 6909055, 44000 6909047)"
@test nonzero_raster_closed_polygon_string(fnas[2]) == "(44001 6909047, 44006 6909047, 44006 6909055, 44001 6909055, 44001 6909047)" 
# Lower level tests
fna = fnas[2]
g = readclose(fna)
z = g.A[:,:,1]
Irng = BitmapMaps.CartesianIndices(BitmapMaps.unpadded_indices(z))
@test sum(Float64.(z[Irng])) == sum(Float64.(z))



cd(olddir)


