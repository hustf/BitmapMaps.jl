using Test
using BitmapMaps
# BM paritioning corresponds to resource/map_sheet_utm_pix.svg
pix_width, pix_height, sheet_pix_width, sheet_pix_height, sw_corner, pix_to_utm_factor = (9, 8, 3, 4, (43999, 6909048), 2)
bmp = BmPartition(pix_width, pix_height, sheet_pix_width, sheet_pix_height, sw_corner, pix_to_utm_factor, "BitmapMaps/test")


for (p, zfi) in zip([bmp[1,1], bmp[2,2]], ["../resource/eksport_796345_20240420.zip", "../resource/eksport_796340_20240420.zip"])
    BitmapMaps.establish_folder(p)
    zipfi = joinpath(@__DIR__, zfi)
    dest = joinpath(BitmapMaps.full_folder_path(p), splitdir(zipfi)[2])
    if ! isfile(dest)
        cp(zipfi, dest)
    end
    # A downloaded and relevant zip file now exists as if downloaded by user.
    # Example export settings are found in /resource/hoydedata export settings.png
    # Other export options work, too.
    BitmapMaps.unzip_tif(p)
    if p == bmp[2,2]
        @test isfile(joinpath(BitmapMaps.full_folder_path(p), "dom1", "data", "dom1-33-1-428-189-63.tif"))
    end
    if p == bmp[1, 1]
        @test isfile(joinpath(BitmapMaps.full_folder_path(p), "dom1", "data", "dom1-33-1-428-189-53.tif"))
        @test isfile(joinpath(BitmapMaps.full_folder_path(p), "dom1", "data", "dom1-33-1-428-189-63.tif"))
    end
end




# Now consolidate
for p in [bmp[1,1], bmp[2,2]]
    BitmapMaps.consolidate_elevation_data(p)
    @test isfile(joinpath(BitmapMaps.full_folder_path(p), BitmapMaps.CONSOLIDATED_FNAM))
end

BitmapMaps.consolidate_elevation_data(bmp)