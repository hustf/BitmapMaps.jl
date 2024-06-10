# The algorithm and parameter values is optimized for in environments/tutorial_images/image_segmentation.jl
# Since these parameter values are well optimized, we use those hardcoded here.
#
# Work with the limited data in /resource. 

using Test
using BitmapMaps

###################################################
# Preparation (much the same as in `test_pipeline`)
###################################################

pth = "BitmapMaps\\test4"
# cleanup 
if ispath(joinpath(homedir(), pth))
    sleep(1) # prevent resource busy error....
    rm(joinpath(homedir(), pth), recursive = true)
end

# Work with limited data from /resource
# Unzip the files to a temporary folder, where the folder name does not provide relevant info. 
tmpdir_topo_relief = mktempdir()

let
    zfi = "../resource/eksport_800912_20240512.zip"
    zipfi = joinpath(@__DIR__, zfi)
    dest = joinpath(tmpdir_topo_relief, splitdir(zipfi)[2])
    if ! isfile(dest)
        cp(zipfi, dest)
    end
    println(dest)
end


# A zip file now exists in tmpdir_topo_relief, as if downloaded by user.
# Extract and inspect. 
unzip_tif(tmpdir_topo_relief)
fna = first(tif_full_filenames_buried_in_folder(tmpdir_topo_relief))

@test nonzero_raster_closed_polygon_string(fna)== "POLYGON ((18294 6937562, 18449 6937562, 18449 6937717, 18294 6937717, 18294 6937562))"


##################
# More preparation
##################

smb = SheetMatrixBuilder(   155, # bm_width_cell
                            155, # bm_height_cell
                             155, # sheet_width_cell
                             155, # sheet_height_cell
                               1, # nrows
                               1, # ncols
                (18294, 6937562), # southwest_corner
    CartesianIndices((1:1, 1:1)), # sheet_indices
                               1, # cell_to_utm_factor
                             189, # sheet_width_mm
                             pth) # pth

@test size(smb[1].cell_iter) == ( sheet_height_cell(smb), sheet_width_cell(smb))

# Establish test folder hierarchy
@test BitmapMaps.establish_folder.(smb) == [true]

# Now copy the relevant file to the relevant directory
fnas = copy_relevant_tifs_to_folder(tmpdir_topo_relief, smb)
# Consolidate
@test BitmapMaps.consolidate_elevation_data.(smb) == [true]
for sb in smb
    fna = joinpath(full_folder_path(sb), BitmapMaps.CONSOLIDATED_FNAM)
    # Get elevation matrix
    za = let 
        z = readclose(fna)
        transpose(z.A[:, :, 1])
    end
    @test sum(za) > 0
end


##############################
# Start working on topo_relief
##############################
sb = smb[1]
BitmapMaps.topo_relief(sb)

outfnam = joinpath(full_folder_path(sb), BitmapMaps.TOPORELIEF_FNAM)
isfile(outfnam)
rm(outfnam)
BitmapMaps.topo_relief(sb)




sb.cell_iter









SheetMatrixBuilder(;
        bm_width_cell                     = 1592,
        bm_height_cell                    = 1150,
        sheet_width_cell                  = 796,
        sheet_height_cell                 = 1150,
        nrows                            = 1,
        ncols                            = 2,
        southwest_corner                 = (43200, 6909000),
        sheet_indices                    = CartesianIndices((1:1, 1:2)),
        cell_to_utm_factor                = 1,
        pth                              = "BitmapMaps\test3")

SheetMatrixBuilder(
        1592,
        1150,
        796,
        1150,
        1,
        2,
        (43200, 6909000),
        CartesianIndices((1:1, 1:2)),
         1,
        "BitmapMaps\test3")












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
unzip_tif(tmpdir_pipeline)
fnas = tif_full_filenames_buried_in_folder(tmpdir_pipeline)
@test nonzero_raster_closed_polygon_string.(fnas)== ["POLYGON ((43200 6909000, 44000 6909000, 44000 6909600, 43200 6909600, 43200 6909000))", "POLYGON ((44000 6909000, 44800 6909000, 44800 6909600, 44000 6909600, 44000 6909000))"]

# Calculate widths, densities and so forth to make a bitmapmap which uses the entire width of the data we have here and fills two A4 sheets
# Since we're dealing with integers, it's a bit of an iteration to get it right
pth = "BitmapMaps\\test3"
nrc = (1, 2)
cell_to_utm_factor = 1
data_cell_width = 44800 - 43200
sheet_width_cell = Int(round(data_cell_width / 2)) 
pwi_max_inch = 191 / 25.4 # mm / (mm / inch) - sheet width in inches
density_pt_m⁻¹ = Int(ceil(sheet_width_cell / pwi_max_inch ))  # Pixels per inch so as to fit roughly sheet_width_cell pixels or more on a sheet
# Reduce the printable width a little from the maximum, so as to fit close to sheet_width_cell pixels
pwi_inch = sheet_width_cell / density_pt_m⁻¹ # cells / (cells / inch)
sheet_width_mm =  Int(floor(pwi_inch * 25.4))  # mm = inch  * (mm/ inch)
sheet_height_mm = Int(ceil((275 / 191) * sheet_width_mm))
complete_sheets_first = false
# Let the pipeline establish the folder structure first:
smb = @test_logs(
    (:info, r"No .tif files"),
    (:info, r"Could not make consolidated"),
    (:warn, r"Could not finish consolidate_elevation_data"),
    run_bitmapmap_pipeline(;nrc, cell_to_utm_factor, pth, southwest_corner = (43200, 6909000), 
    density_pt_m⁻¹, sheet_width_mm, sheet_height_mm, complete_sheets_first)
    )



@test abs(sheet_width_cell(smb) - sheet_width_cell) <= 4 # Good enough, fits inside data
@test nrows(smb) == 1
@test BitmapMaps.ncols(smb) == 2
@test ispath(full_folder_path(smb[1,1]))
@test ispath(full_folder_path(smb[1,2]))

# Now copy the relevant files to the relevant directories
fnas = copy_relevant_tifs_to_folder(tmpdir_pipeline, smb)
# One of the sheets is based off two files, the other sheet is fully covered by one of those two files.
@test length(fnas) == 3

# Run the pipeline again, this time proceeding further
smb = run_bitmapmap_pipeline(;nrc, cell_to_utm_factor, pth, southwest_corner = (43200, 6909000), 
    density_pt_m⁻¹, sheet_width_mm, sheet_height_mm)
@test ispath(joinpath(full_folder_path(smb[1,1]), BitmapMaps.CONSOLIDATED_FNAM))
@test ispath(joinpath(full_folder_path(smb[1,2]), BitmapMaps.CONSOLIDATED_FNAM))


# cleanup 
if ispath(joinpath(homedir(), pth))
    sleep(1) # prevent resource busy error....
    rm(joinpath(homedir(), pth), recursive = true)
end

