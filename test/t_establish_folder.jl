using Test
using BitmapMaps

for fo in ["test", "te st", "1\\2", "3 4 / 5", "..", "a/../b/c/d", ""]
    # BM paritioning corresponds to resource/map_sheet_utm_pix.svg
    pix_width, pix_height, sheet_pix_width, sheet_pix_height, sw_corner, pix_to_utm_factor = (9, 8, 3, 4, (43999, 6909048), 2)
    bmp = BmPartition(pix_width, pix_height, sheet_pix_width, sheet_pix_height, sw_corner, pix_to_utm_factor, fo)
    p = bmp[1]
    fofo = joinpath(homedir(), p.pthsh)
    @test (1, 1) == BitmapMaps.parse_folder_name(fofo)[1:2]
    @test southwest_corner(p) == BitmapMaps.parse_folder_name(fofo)[3:4]
    @test northeast_external_corner(p) == BitmapMaps.parse_folder_name(fofo)[5:6]
end
