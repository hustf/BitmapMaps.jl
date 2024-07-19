using Test
using BitmapMaps

for fo in ["test", "te st", "1\\2", "3 4 / 5", "..", "a/../b/c/d", ""]
    # SheetMatrixBuilder corresponds to resource/matrix_sheet_cell_utm.svg
    width_cell, height_cell, sh_width_cell, sheet_height_cell, sw_corner, cell_to_utm_factor = (9, 8, 3, 4, (44000, 6909047), 2)
    smb = SheetMatrixBuilder(width_cell, height_cell, sh_width_cell, sheet_height_cell, sw_corner, cell_to_utm_factor, fo)
    p = smb[1]
    fofo = joinpath(homedir(), p.pthsh)
    @test (1, 1) == BitmapMaps.parse_folder_name(fofo)[1:2]
    @test southwest_corner(p) == BitmapMaps.parse_folder_name(fofo)[3:4]
    @test northeast_external_corner(p) == BitmapMaps.parse_folder_name(fofo)[5:6]
end

