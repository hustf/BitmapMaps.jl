using Test
using BitmapMaps
# SheetBuilder corresponds to resource/matrix_sheet_pix_utm.svg
bm_pix_width, bm_pix_height, sheet_pix_width, sheet_pix_height, sw_corner, pix_to_utm_factor = (9, 8, 3, 4, (43999, 6909048), 2)
smb = SheetMatrixBuilder(bm_pix_width, bm_pix_height, sheet_pix_width, sheet_pix_height, sw_corner, pix_to_utm_factor, "BitmapMaps/test")



for (p, zfi) in zip([smb[1,1], smb[2,2]], ["../resource/eksport_796345_20240420.zip", "../resource/eksport_796340_20240420.zip"])
    BitmapMaps.establish_folder(p)
    zipfi = joinpath(@__DIR__, zfi)
    dest = joinpath(BitmapMaps.full_folder_path(p), splitdir(zipfi)[2])
    if ! isfile(dest)
        cp(zipfi, dest)
    end
    # A downloaded and relevant zip file now exists as if downloaded by user.
    # Example export settings are found in /resource/hoydedata export settings.png
    # Other export options work, too.
    unzip_tif(p)
    if p == smb[2,2]
        @test isfile(joinpath(BitmapMaps.full_folder_path(p), "dom1", "data", "dom1-33-1-428-189-63.tif"))
    end
    if p == smb[1, 1]
        @test isfile(joinpath(BitmapMaps.full_folder_path(p), "dom1", "data", "dom1-33-1-428-189-53.tif"))
        @test isfile(joinpath(BitmapMaps.full_folder_path(p), "dom1", "data", "dom1-33-1-428-189-63.tif"))
    end
end




# Now consolidate (multiple) source files into CONSOLIDATED_FNAM in each sheet's folder.
for p in [smb[2,2], smb[1,1] ]
    fnam_out = joinpath(BitmapMaps.full_folder_path(p), BitmapMaps.CONSOLIDATED_FNAM)
    if isfile(fnam_out)
        rm(fnam_out)
    end
    BitmapMaps.consolidate_elevation_data(p)
    @test isfile(fnam_out)
end

# Clean up
for p in [smb[2,2], smb[1,1] ]
    fnam_out = joinpath(BitmapMaps.full_folder_path(p), BitmapMaps.CONSOLIDATED_FNAM)
    if isfile(fnam_out)
        sleep(1) # prevent resource busy error....
        rm(fnam_out)
    end
end

