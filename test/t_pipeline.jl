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
    run_bitmapmap_pipeline(;nrc = (1, 1), pix_to_utm_factor = 1, pth = "BitmapMaps\\test2"))

@test ispath(joinpath(homedir(), "BitmapMaps", "test2"))




# Work with the limited data in /resource
# Unzip the files to a temporary folder, where the folder name does not provide relevant info. 
tmpdir_pipeline = mktempdir()
let
    for zfi in ["../resource/eksport_796345_20240420.zip", "../resource/eksport_796340_20240420.zip"]
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
unzip_tif(tmpdir_pipeline)
fnas = tif_full_filenames_buried_in_folder(tmpdir_pipeline)
@test BitmapMaps.bounding_box_polygon_string.(fnas)== ["POLYGON ((43200 6909000, 44000 6909000, 44000 6909600, 43200 6909600, 43200 6909000))", "POLYGON ((44000 6909000, 44800 6909000, 44800 6909600, 44000 6909600, 44000 6909000))"]

# Calculate widths, densities and so forth to make a bitmapmap which uses the entire width of the data we have here and fills two A4 sheets
# Since we're dealing with integers, it's a bit of an iteration to get it right
pth = "BitmapMaps\\test3"
nrc = (1, 2)
pix_to_utm_factor = 1
data_pix_width = 44800 - 43200
sheet_pix_width = Int(round(data_pix_width / 2)) 
pwi_max_inch = 191 / 25.4 # mm / (mm / inch) - sheet width in inches
pdens_dpi = Int(ceil(sheet_pix_width / pwi_max_inch ))  # Pixels per inch so as to fit roughly sheet_pix_width pixels or more on a sheet
# Reduce the printable width a little from the maximum, so as to fit close to sheet_pix_width pixels
pwi_inch = sheet_pix_width / pdens_dpi # pix / (pix / inch)
pwi =  Int(floor(pwi_inch * 25.4))  # mm = inch  * (mm/ inch)
phe = Int(ceil((275 / 191) * pwi))
complete_sheets_first = false
# Let the pipeline establish the folder structure first:
smb = @test_logs(
    (:info, r"No .tif files"),
    (:info, r"Could not make consolidated"),
    (:warn, r"Could not finish consolidate_elevation_data"),
    run_bitmapmap_pipeline(;nrc, pix_to_utm_factor, pth, southwest_corner = (43200, 6909000), 
    pdens_dpi, pwi, phe, complete_sheets_first)
    )



@test abs(smb.sheet_pix_width - sheet_pix_width) <= 4 # Good enough, fits inside data
@test smb.nrows == 1
@test smb.ncols == 2
@test ispath(full_folder_path(smb[1,1]))
@test ispath(full_folder_path(smb[1,2]))

# Now copy the relevant files to the relevant directories
fnas = copy_relevant_tifs_to_folder(tmpdir_pipeline, smb)
# One of the sheets is based off two files, the other sheet is fully covered by one of those two files.
@test length(fnas) == 3

# Run the pipeline again, this time proceeding further
smb = run_bitmapmap_pipeline(;nrc, pix_to_utm_factor, pth, southwest_corner = (43200, 6909000), 
    pdens_dpi, pwi, phe)
@test ispath(joinpath(full_folder_path(smb[1,1]), BitmapMaps.CONSOLIDATED_FNAM))
@test ispath(joinpath(full_folder_path(smb[1,2]), BitmapMaps.CONSOLIDATED_FNAM))

# Cleanup
for fo in ["test1", "test2", "test3"]
    if ispath(joinpath(homedir(), "BitmapMaps", fo))
        # reading a .tif files unfortunately leaves them open. 
        # The files are released when we close Julia, so not worth putting effort into fixing. 
        # We can cleanup at start instead.
        try
            rm(joinpath(homedir(), "BitmapMaps", fo), recursive = true)
        catch
        end
    end
end
