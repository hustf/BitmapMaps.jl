module BitmapMaps
# Extend methods for builders
import Base: iterate, length, getindex, size, axes
# Read defaults
using IniFile
# For loading and displaying png:
import FileIO
import ImageShow
using FileIO: load
# For adding physical print widths to png files in png_phys.jl
import PNGFiles
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
using LuxorLayout: Point # Maybe no need, but we won't discard compatibility easily.
import Base: show
# ZipFile and comparing file hashes:
import ZipFile
import SHA
using SHA: sha256
# GeoArray, hypsometric colors
import GeoArrays
using GeoArrays: GeoArray, bbox, bbox_overlap, Vertex, coords, indices, crop
# Feedback
import Dates
# Calculating gradients, topo map
import ImageFiltering
using ImageFiltering: imgradients, KernelFactors, mapwindow, mapwindow!, Kernel, imfilter, centered
# Pick color from lightness
import Colors
import ColorTypes
using ColorTypes: RGBA, RGB, AbstractGray, AbstractRGB, Gray, red, green, blue # Todo clean
# Identify water
import ImageSegmentation
using ImageSegmentation: felzenszwalb, labels_map, segment_pixel_count, segment_mean
import ImageMorphology
using ImageMorphology: erode!, dilate!, dilate, GuoAlgo, thinning
# Grid
import Geodesy
using Geodesy: LLAfromUTMZ, UTMZfromLLA, wgs84, UTMZ, UTM
# Elevation contours, consolidate feedback
using ColorTypes: GrayA, gray, alpha
import ImageCore
using ImageCore: channelview, scaleminmax, colorview, N0f8
using ImageFiltering: FIRTiled, Fill
# Summit prominence
using ImageMorphology: MaxTree
import DelimitedFiles
using DelimitedFiles: readdlm, writedlm
# Overlays to composite
import ColorBlendModes
using ColorBlendModes: CompositeDestinationOver, BlendLighten
# Exports
export save_png_with_phys
export run_bitmapmap_pipeline
# Export of builders and and utilties
export SheetBuilder, SheetMatrixBuilder, full_folder_path
export northeast_corner, northeast_external_corner, northeast_internal_corner
export northwest_corner
export southeast_corner, southeast_external_corner, southeast_internal_corner
export southwest_corner, southwest_external_corner, southwest_internal_corner

export geo_grid_centre_single, geo_centre, bbox_external_string, polygon_string, show_augmented
export geo_area, nonzero_raster_closed_polygon_string, raster_polygon_string, nonzero_raster_rect
export copy_relevant_tifs_to_folder, tif_full_filenames_buried_in_folder, unzip_tif
export show_derived_properties, readclose, display_if_vscode
#
# Not exported because it's too generic: mark_at! and line!

const CONSOLIDATED_FNAM = "Consolidated.tif"
const TOPORELIEF_FNAM = "Toporelief.png"
const WATER_FNAM = "Water.png"
const CONTOUR_FNAM = "Contour.png"
const GRID_FNAM = "Grid.png"
const RIDGE_FNAM = "Ridge.png"
const MARKERS_FNAM = "Markers.png"
const MAX_ELEVATION_ABOVE_FNAM = "_Max_elevation_above.mat"
const SUMMITS_FNAM = "Summits.csv"
const COMPOSITE_FNAM = "Composite.png"

include("ini_file.jl")
include("png_phys.jl")
include("filter.jl")
include("pallette.jl")
include("builders.jl")
include("builders_utilties.jl")
include("geoarray_utilties.jl")
include("topo_relief.jl")
include("pipeline.jl")
include("establish_folder.jl")
include("unzip_tif.jl")
include("consolidate_elevation_data.jl")
include("user_utilties.jl")
include("water.jl")
include("contour.jl")
include("grid.jl")
include("ridges.jl")
include("layers.jl")
include("mark_utils.jl")
include("summit_markers.jl")

end