using Test
using BitmapMaps


# Corresponds to resource/matrix_sheet_cell_utm.svg
sw_corner, nrc, cell2utm = (44000, 6909047), (2,3), 2
sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, pth = 3, 4, 1000, "BitmapMaps/test_unzip"
smb = SheetMatrixBuilder(sw_corner, nrc, cell2utm, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, pth)

# Pre-test clean up
rm(full_folder_path(smb), force = true, recursive = true)

# Copy .tif files into homedir() / BitmapMaps / test_unzip.
for (p, zfi) in zip([smb[1,1], smb[2,2]], ["../resource/eksport_796345_20240420.zip", "../resource/eksport_796340_20240420.zip"])
    BitmapMaps.establish_folder(p)
    zipfi = abspath(joinpath(@__DIR__, zfi))
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
        @test polygon_string(fn1) == "MULTIPOLYGON (\n                   ((43200 6909000, 44000 6909000, 44000 6909600, 43200 6909600, 43200 6909000)),\n                   ((43999 6909048, 44000 6909048, 44000 6909056, 43999 6909056, 43999 6909048)))"
        @test polygon_string(fn2) == "MULTIPOLYGON (\n                   ((44000 6909000, 44800 6909000, 44800 6909600, 44000 6909600, 44000 6909000)),\n                   ((44000 6909048, 44005 6909048, 44005 6909056, 44000 6909056, 44000 6909048)))"
    else
        fn = joinpath(full_folder_path(p), "dom1", "data", "dom1-33-1-428-189-63.tif") # Same file name as fn2, not same content.
        @test isfile(fn)
        @test polygon_string(fn) == "MULTIPOLYGON (\n                   ((44000 6909000, 44800 6909000, 44800 6909600, 44000 6909600, 44000 6909000)),\n                   ((44005 6909056, 44011 6909056, 44011 6909064, 44005 6909064, 44005 6909056)))"
    end
end


# Test constructing a geoarrray corresponding to a sheet builder.
@testset "A" begin
    p = smb[2, 2]
    fofo = full_folder_path(smb[2, 2])
    r, c, min_x, min_y, max_x, max_y = BitmapMaps.parse_folder_name(fofo)
    w = max_x - min_x
    h = max_y - min_y
    g = let
        A = zeros(Float32, w, h)
        linear = BitmapMaps.SA[1.0 0.0; 0.0 -1.0]
        translation = BitmapMaps.SA[Float64(min_x), Float64(max_y)]
        f = BitmapMaps.AffineMap(linear, translation)
        BitmapMaps.GeoArray(A, f)
    end
    fill!(g.A, Float32(1))
    # This expression makes two boxes: one with all data, and one with non-zero data.
    # Here, both boxes are similar, but the check for that does not work, because we
    # did not fully make the transition to using the new module Extents here.
    @test polygon_string(g) == "MULTIPOLYGON (\n                   ((44006 6909055, 44012 6909055, 44012 6909063, 44006 6909063, 44006 6909055)),\n                   ((44006 6909055, 44012 6909055, 44012 6909063, 44006 6909063, 44006 6909055)))"
    @test polygon_string(smb[2, 2]) == "MULTIPOLYGON (\n                   ((44006 6909055, 44012 6909055, 44012 6909063, 44006 6909063, 44006 6909055)))"
    g.A[:, 1, 1] .= Float32(0)
    #     
    @test polygon_string(g) == "MULTIPOLYGON (\n                   ((44006 6909055, 44012 6909055, 44012 6909063, 44006 6909063, 44006 6909055)),\n                   ((44006 6909055, 44012 6909055, 44012 6909062, 44006 6909062, 44006 6909055)))"
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
    @test sum(g.A[:, :]) > 0
    # We lack elevation data for two sides
    @test sum(iszero.(g)) == 13
    if p == smb[2, 2]
        @test polygon_string(g) == "MULTIPOLYGON (\n                   ((44006 6909055, 44012 6909055, 44012 6909063, 44006 6909063, 44006 6909055)),\n                   ((44006 6909056, 44011 6909056, 44011 6909063, 44006 6909063, 44006 6909056)))"
    else
        @test polygon_string(g) == "MULTIPOLYGON (\n                   ((44000 6909047, 44006 6909047, 44006 6909055, 44000 6909055, 44000 6909047)),\n                   ((44000 6909048, 44005 6909048, 44005 6909055, 44000 6909055, 44000 6909048)))"
    end
end

# Clean up
rm(full_folder_path(smb), force = true, recursive = true)


# Translate the builder so that files cover the two sheets tested.
sw_corner, nrc, cell2utm = (43999, 6909048)  , (2,3), 2
sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, pth = 3, 4, 1000, "BitmapMaps/test_unzip"
smb = SheetMatrixBuilder(sw_corner, nrc, cell2utm, sheet_width_mm, sheet_height_mm, density_pt_m⁻¹, pth)

# Copy .tif files into homedir() / BitmapMaps / test_unzip.
# The exact same files are copied this time.
for (p, zfi) in zip([smb[1,1], smb[2,2]], ["../resource/eksport_796345_20240420.zip", "../resource/eksport_796340_20240420.zip"])
    BitmapMaps.establish_folder(p)
    zipfi = abspath(joinpath(@__DIR__, zfi))
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
        @test polygon_string(fn1) == "MULTIPOLYGON (\n                   ((43200 6909000, 44000 6909000, 44000 6909600, 43200 6909600, 43200 6909000)),\n                   ((43999 6909048, 44000 6909048, 44000 6909056, 43999 6909056, 43999 6909048)))"
        @test polygon_string(fn2) == "MULTIPOLYGON (\n                   ((44000 6909000, 44800 6909000, 44800 6909600, 44000 6909600, 44000 6909000)),\n                   ((44000 6909048, 44005 6909048, 44005 6909056, 44000 6909056, 44000 6909048)))"
    else
        fn = joinpath(full_folder_path(p), "dom1", "data", "dom1-33-1-428-189-63.tif") # Same file name as fn2, not same content.
        @test isfile(fn)
        @test polygon_string(fn) == "MULTIPOLYGON (\n                   ((44000 6909000, 44800 6909000, 44800 6909600, 44000 6909600, 44000 6909000)),\n                   ((44005 6909056, 44011 6909056, 44011 6909064, 44005 6909064, 44005 6909056)))"
    end
end

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
    @test sum(g.A[:, :]) > 0
    # We have full elevation data
    @test sum(iszero.(g)) == 0
    # Broken due to incomplete implementation of Extents.Extent.
    @test_broken polygon_string(g) == polygon_string(p)
    # 
    @test southwest_external_corner(g) == southwest_external_corner(p)
    @test northeast_external_corner(g) == northeast_external_corner(p)

end

# Cleanup
if ispath(joinpath(homedir(), pth))
    rm(joinpath(homedir(), pth), recursive = true)
end
