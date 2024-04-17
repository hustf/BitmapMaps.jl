using Test
using BitmapMaps

run_bitmapmap_pipeline(;pth = "BitmapMaps\\test1")
run_bitmapmap_pipeline(;nrc = (1, 1), pix_to_utm_factor = 1, pth = "BitmapMaps\\test2");
@test true