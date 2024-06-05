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
withenv("JULIA_DEBUG" => "BitmapMaps") do
    @test_logs(
        (:debug, r"Name and path similarity, made unique file namme."),
        unzip_tif(tmpdir_pipeline))
end



fnas = tif_full_filenames_buried_in_folder(tmpdir_pipeline)
@test polygon_string.(fnas)== ["POLYGON ((44000 6909047, 44001 6909047, 44001 6909055, 44000 6909055, 44000 6909047))", "POLYGON ((44001 6909047, 44006 6909047, 44006 6909055, 44001 6909055, 44001 6909047))", "POLYGON ((44006 6909055, 44012 6909055, 44012 6909063, 44006 6909063, 44006 6909055))"]

# Calculate widths, densities and so forth to make a bitmapmap which uses the entire width of the data we have here and fills two A4 sheets
# Since we're dealing with integers, it's a bit of an iteration to get it right
pth = "BitmapMaps\\test3"
nrc = (1, 2)
cell_to_utm_factor = 1
southwest_corner = (44000, 6909047)
data_cell_width = 44012 - southwest_corner[1]
sheet_cell_width = Int(round(data_cell_width / 2)) 
pwi_max_inch = 191 / 25.4 # mm / (mm / inch) - sheet width in inches
pdens_dpi = Int(ceil(sheet_cell_width / pwi_max_inch ))  # Pixels per inch so as to fit roughly sheet_cell_width pixels or more on a sheet
# Reduce the printable width a little from the maximum, so as to fit close to sheet_cell_width pixels
pwi_inch = sheet_cell_width / pdens_dpi # cells / (cells / inch)
pwi =  Int(floor(pwi_inch * 25.4))  # mm = inch  * (mm/ inch)
phe = Int(ceil((275 / 191) * pwi))  # Height from A4 aspect ratio
complete_sheets_first = false
# Let the pipeline establish the folder structure first. Data files are missing, so will fail gracefully.
smb = @test_logs(
    (:info, r"No .tif files"),
    (:info, r"Could not make consolidated"),
    (:warn, r"Could not finish consolidate_elevation_data"),
    run_bitmapmap_pipeline(;nrc, cell_to_utm_factor, pth, southwest_corner, 
    pdens_dpi, pwi, phe, complete_sheets_first)
    )
@test abs(smb.sheet_cell_width - sheet_cell_width) <= 4 # Good enough, fits inside data
@test smb.nrows == 1
@test smb.ncols == 2
@test ispath(full_folder_path(smb[1,1]))
@test ispath(full_folder_path(smb[1,2]))

# Now copy the relevant, unzipped files to the directories of each sheet
fnas = copy_relevant_tifs_to_folder(tmpdir_pipeline, smb)
# Both sheets have two input files. One file is duplicated to both sheet's folders.
@test length(fnas) == 4

# Run the pipeline again, this time proceeding further
# TODO: Why does smb1 differ from smb? sheet_cell_width changes! New folders are made!
smb1 = run_bitmapmap_pipeline(smb)
@test smb1 == smb
smb = run_bitmapmap_pipeline(;nrc, cell_to_utm_factor, pth, southwest_corner, 
    pdens_dpi, pwi, phe, )

@test ispath(joinpath(full_folder_path(smb[1,1]), BitmapMaps.CONSOLIDATED_FNAM))
@test ispath(joinpath(full_folder_path(smb[1,2]), BitmapMaps.CONSOLIDATED_FNAM))

# Cleanup
for fo in ["test1", "test2", "test3"]
    if ispath(joinpath(homedir(), "BitmapMaps", fo))
        rm(joinpath(homedir(), "BitmapMaps", fo), recursive = true)
    end
end
