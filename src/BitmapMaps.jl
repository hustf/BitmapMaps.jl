module BitmapMaps
using IniFile
import Base: iterate, length, getindex, size, axes
# For loading and displaying png:
import FileIO
import ImageShow
using FileIO: load

# For adding physical print widths to png files in png_phys.jl
import PNGFiles
using PNGFiles: N0f8, RGBA, AbstractGray, AbstractRGB
using PNGFiles: _save, png_init_io, close_png, png_uint_32
using PNGFiles: Z_BEST_SPEED, Z_RLE, PNG_FILTER_PAETH, Z_DEFAULT_STRATEGY
using PNGFiles: Z_FIXED, Z_NO_COMPRESSION, Z_BEST_COMPRESSION
using PNGFiles: create_write_struct, create_info_struct, png_set_pHYs
# For inspecting pHYs chunk in png
using PNGFiles: open_png, create_read_struct, PNG_COLOR_TYPE_RGB
using PNGFiles: png_set_sig_bytes, PNG_BYTES_TO_CHECK, png_read_info, png_get_image_width, png_get_image_height
using PNGFiles: png_get_color_type, png_get_bit_depth, png_color_16p, png_get_bKGD, PNG_COLOR_TYPE_PALETTE
using PNGFiles: PNG_COLOR_TYPE_GRAY, png_destroy_read_struct, close_png
using PNGFiles: PNG_INFO_pHYs, png_get_pHYs
# For sheet builders
import LuxorLayout
using LuxorLayout: Point # Maybe no need
import Base: show
# ZipFile
import ZipFile
# GeoArray
import GeoArrays
using GeoArrays: GeoArray, bbox, bbox_overlap, Vertex, coords, indices, crop
# Feedback
import Dates
# Exports
export save_png_with_phys
export run_bitmapmap_pipeline
# Export of builders and and utilties
export SheetBuilder, SheetMatrixBuilder, row_col_of_sheet
export northeast_corner, northeast_external_corner, northeast_internal_corner
export northwest_corner
export southeast_corner, southeast_external_corner, southeast_internal_corner
export southwest_corner, southwest_external_corner, southwest_internal_corner

export geo_grid_centre_single, geo_centre, bounding_box_external_string, show_augmented
export geo_area
export copy_relevant_tifs_to_folder
export show_augmented_properties

const CONSOLIDATED_FNAM = "Consolidated.tif"

include("ini_file.jl")
include("png_phys.jl")
include("builders.jl")
include("builders_utilties.jl")
include("pipeline.jl")
include("establish_folder.jl")
include("unzip_tif.jl")
include("consolidate_elevation_data.jl")
include("user_utilties.jl")
end