using Test
using BitmapMaps

pwi = BitmapMaps.get_config_value("Printer consistent capability", "Printable width mm", Int)
@test pwi isa Int
@test pwi > 0
