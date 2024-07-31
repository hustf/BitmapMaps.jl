using Test
using BitmapMaps
using BitmapMaps: get_kw_or_config_value
@test BitmapMaps.get_config_value("Printer consistent capability", "Printable width mm", Int) isa Int
@test BitmapMaps.get_config_value("Printer consistent capability", "Printable width mm", Int) > 0
@test BitmapMaps.get_config_value("Geographical area",
    "Output density, i.e. 'cells' / 'dots' / 'points' or 'pixels' per paper meter", Int) isa Int
@test get_kw_or_config_value(:sheet_width_mm ,"Printer consistent capability", "Printable width mm", Int) isa Int
@test get_kw_or_config_value(:sheet_height_mm ,"Printer consistent capability", "Printable height mm", Int) isa Int
# Other types than Int...
@test get_kw_or_config_value(:southwest_corner ,"Geographical area", "Southwest corner (utm easting northing)", Tuple{Int, Int}) isa Tuple{Int, Int}
@test get_kw_or_config_value(:pth, "File folder", "Top folders path under homedir()", String) isa String
