using Test
using BitmapMaps

img = BitmapMaps.load(joinpath(@__DIR__, "..", "resource", "bitmap_detail.png"))
tmpdir_png = mktempdir()
fna = joinpath(tmpdir_png, "600_dpi.png")
save_png_with_phys(fna, img, Int(round(600/0.0254)))
res_x, res_y, unit_type = BitmapMaps.get_pHYs_chunk_res_x_y_unit(fna; silent = true)
@test round(res_x * 0.0254) == 600
# cleanup
if ispath(joinpath(homedir(), tmpdir_png))
    rm(joinpath(homedir(), tmpdir_png), recursive = true)
end
