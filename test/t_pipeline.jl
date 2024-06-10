using Test
using BitmapMaps

# cleanup
for fo in ["test1", "test2", "test3"]
    if ispath(joinpath(homedir(), "BitmapMaps", fo))
        sleep(1) # prevent resource busy error....
        rm(joinpath(homedir(), "BitmapMaps", fo), recursive = true)
    end
end

# First some impossible jobs to finish. 
@test_logs(
    (:info, r"Since geographical area per sheet is  > 16km²"),
    (:info, r"No .tif files"),
    (:info, r"Could not make consolidated"),
    (:warn, r"Could not finish consolidate_elevation_data"),
    run_bitmapmap_pipeline(;pth = "BitmapMaps\\test1"))


@test ispath(joinpath(homedir(), "BitmapMaps", "test1"))

@test_logs(
    (:info, r"No .tif files"),
    (:info, r"Could not make consolidated"),
    (:warn, r"Could not finish consolidate_elevation_data"),
    #match_mode = :any,
    run_bitmapmap_pipeline(;nrc = (1, 1), cell_to_utm_factor = 1, pth = "BitmapMaps\\test2"))

@test ispath(joinpath(homedir(), "BitmapMaps", "test2"))


# Work with the limited data in /resource
# Unzip the files to a temporary folder, where the folder name does not provide relevant info. 
tmpdir_pipeline = mktempdir()
let
    for zfi in ["../resource/eksport_796345_20240420.zip", 
            "../resource/eksport_796340_20240420.zip", 
            "../resource/eksport_826662_20240610.zip"]
        zipfi = joinpath(@__DIR__, zfi)
        dest = joinpath(tmpdir_pipeline, splitdir(zipfi)[2])
        if ! isfile(dest)
            cp(zipfi, dest)
        end
        println(dest)
    end
end
# A zip file now exists in tmpdir_pipeline, as if downloaded by user.
# Extract and inspect. 
withenv("JULIA_DEBUG" => "BitmapMaps") do
    @test_logs(
        (:debug, r"Name and path similarity, made unique file namme."),
        (:debug, r"Name and path similarity, made unique file namme."),
        unzip_tif(tmpdir_pipeline))
end



fnas = tif_full_filenames_buried_in_folder(tmpdir_pipeline)
@test polygon_string.(fnas) == ["POLYGON ((44000 6909047, 44001 6909047, 44001 6909055, 44000 6909055, 44000 6909047))", 
    "POLYGON ((44001 6909047, 44006 6909047, 44006 6909055, 44001 6909055, 44001 6909047))", 
    "POLYGON ((44007 6909046, 44013 6909046, 44013 6909054, 44007 6909054, 44007 6909046))", 
    "POLYGON ((44006 6909055, 44012 6909055, 44012 6909063, 44006 6909063, 44006 6909055))"]
nonzero_raster_rect.(fnas) == [(min_x = 44000, min_y = 6909047, max_x = 44001, max_y = 6909055), 
    (min_x = 44001, min_y = 6909047, max_x = 44006, max_y = 6909055), 
    (min_x = 44007, min_y = 6909046, max_x = 44013, max_y = 6909054), 
    (min_x = 44006, min_y = 6909055, max_x = 44012, max_y = 6909063)]
# Calculate widths, densities and so forth to make a bitmapmap which uses the entire width of the data we have here and 
# fills two A4 sheets.
pth = "BitmapMaps\\test3"
nrc = (1, 2)
cell_to_utm_factor = 1
southwest_c = (44000, 6909047)
data_cell_width = 44012 - southwest_c[1]
sheet_width_cell = Int(round(data_cell_width / 2)) 
sheet_width_mm = 191 
density_pt_m⁻¹ = Int(ceil(1000 *sheet_width_cell / sheet_width_mm )) 
# Let the pipeline establish the folder structure first. Data files are missing, so will fail gracefully.
smb = @test_logs(
    (:info, r"No .tif files"),
    (:info, r"Could not make consolidated"),
    (:warn, r"Could not finish consolidate_elevation_data"),
    run_bitmapmap_pipeline(;nrc, cell_to_utm_factor, pth, southwest_corner = southwest_c, density_pt_m⁻¹, sheet_width_mm,
       complete_sheets_first = false))
@test BitmapMaps.sheet_width_cell(smb) == sheet_width_cell 
@test BitmapMaps.nrows(smb) == 1
@test BitmapMaps.ncols(smb) == 2
@test ispath(full_folder_path(smb[1,1]))
@test ispath(full_folder_path(smb[1,2]))

# Now copy the relevant, unzipped files to the directories of each sheet
fnas = copy_relevant_tifs_to_folder(tmpdir_pipeline, smb)
@test length(fnas) == 3

# Run the pipeline again, this time proceeding further
smb1 = run_bitmapmap_pipeline(smb, complete_sheets_first = true)
# What's wrong in the second folder????
@test smb1 == smb
@test ispath(joinpath(full_folder_path(smb[1,1]), BitmapMaps.CONSOLIDATED_FNAM))
@test ispath(joinpath(full_folder_path(smb[1,2]), BitmapMaps.CONSOLIDATED_FNAM))

# Cleanup
for fo in ["test1", "test2", "test3"]
    if ispath(joinpath(homedir(), "BitmapMaps", fo))
        rm(joinpath(homedir(), "BitmapMaps", fo), recursive = true)
    end
end
