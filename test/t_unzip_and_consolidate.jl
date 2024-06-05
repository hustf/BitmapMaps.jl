using Test
using BitmapMaps
# SheetBuilder corresponds to resource/matrix_sheet_cell_utm.svg
bm_cell_width, bm_cell_height, sheet_cell_width, sheet_cell_height, sw_corner, cell_to_utm_factor = (9, 8, 3, 4, (44000, 6909047), 2)
smb = SheetMatrixBuilder(bm_cell_width, bm_cell_height, sheet_cell_width, sheet_cell_height, sw_corner, cell_to_utm_factor, 191, "BitmapMaps/test")


# Copy .tif files into homedir() / BitmapMaps / test.
for (p, zfi) in zip([smb[1,1], smb[2,2]], ["../resource/eksport_796345_20240420.zip", "../resource/eksport_796340_20240420.zip"])
    BitmapMaps.establish_folder(p)
    zipfi = joinpath(@__DIR__, zfi)
    dest = joinpath(full_folder_path(p), splitdir(zipfi)[2])
    if ! isfile(dest)
        cp(zipfi, dest)
    end
    # A downloaded and relevant zip file now exists as if downloaded by user.
    # Example export settings are found in /resource/hoydedata export settings.png
    # Other export options work, too.
    unzip_tif(p)
    if p == smb[1, 1]
        fn1 = joinpath(full_folder_path(p), "dom1", "data", "dom1-33-1-428-189-53.tif")
        @test isfile(fn1)
        fn2 = joinpath(full_folder_path(p), "dom1", "data", "dom1-33-1-428-189-63.tif")
        @test isfile(fn2)
        @test polygon_string(fn1) == "POLYGON ((44000 6909047, 44001 6909047, 44001 6909055, 44000 6909055, 44000 6909047))"
        @test polygon_string(fn2) == "POLYGON ((44001 6909047, 44006 6909047, 44006 6909055, 44001 6909055, 44001 6909047))"
    else
        fn = joinpath(full_folder_path(p), "dom1", "data", "dom1-33-1-428-189-63.tif") # Same file name as fn2, not same content.
        @test isfile(fn)
        @test polygon_string(fn) == "POLYGON ((44006 6909055, 44012 6909055, 44012 6909063, 44006 6909063, 44006 6909055))"
    end
end



# Now consolidate (multiple) source files into CONSOLIDATED_FNAM in each sheet's folder.
for p in [smb[2,2], smb[1,1] ]
    fnam_out = joinpath(full_folder_path(p), BitmapMaps.CONSOLIDATED_FNAM)
    # Just in case clean-up didn't work (was a problem prior to `readclose`. Ref issue GeoArrays.jl 159)
    if isfile(fnam_out)
        rm(fnam_out)
    end
    BitmapMaps.consolidate_elevation_data(p)
    @test isfile(fnam_out)
    # Get consolidated geoarray
    g = readclose(fnam_out)
    # Not empty
    @test sum(g.A[:, :, 1]) > 0
    # We found elevation data for all cells
    @test sum(iszero.(g)) == 0
    @test polygon_string(fnam_out) == polygon_string(p)
end

# Clean up
for p in [smb[2,2], smb[1,1] ]
    fnam_out = joinpath(full_folder_path(p), BitmapMaps.CONSOLIDATED_FNAM)
    if isfile(fnam_out)
        rm(fnam_out)
    end
end

