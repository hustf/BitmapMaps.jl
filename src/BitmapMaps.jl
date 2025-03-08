module BitmapMaps
# Extend methods for builders
import Base: iterate, length, getindex, size, axes
# Read defaults
using IniFile
# For loading and displaying png:
import FileIO
using FileIO: load
import ImageShow
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
import Base: show
# ZipFile and comparing file hashes:
import ZipFile
import SHA
import Serialization
using Serialization: serialize, deserialize
using SHA: sha256
# GeoArrays
import GeoArrays
using GeoArrays: GeoArray, AffineMap, bbox
import Extents
using Extents: bounds
import StaticArrays
using StaticArrays: SMatrix, SVector, SA
# Feedback
import Dates
using Dates: now, format
# Calculating gradients, topo map, artifact detection
import ImageFiltering
using ImageFiltering: imgradients, KernelFactors, mapwindow, mapwindow!, Kernel, imfilter, centered
# Pick color from lightness
import Colors
import ColorTypes
using ColorTypes: RGBA, RGB, XYZ, XYZA, AbstractGray, AbstractRGB, Gray, red, green, blue
# Identify water
import ImageSegmentation
using ImageSegmentation: felzenszwalb, labels_map, segment_pixel_count, segment_mean, SegmentedImage
using ImageSegmentation: fast_scanning, prune_segments, segment_labels
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
# Summit prominence, names, and graph_utilties
using ImageMorphology: MaxTree, strel
import DelimitedFiles
using DelimitedFiles: readdlm, writedlm
import Stadnamn
using Stadnamn: point_names
import MetaGraphsNext
using MetaGraphsNext: loadgraph, savegraph, MetaGraph, DiGraph, labels, nv, ne, add_edge!, MGFormat, SimpleDiGraph
using MetaGraphsNext: all_neighbors, neighbor_labels, inneighbor_labels, outneighbor_labels, all_neighbor_labels
using MetaGraphsNext: neighbors, indegree, outdegree, rem_edge!, rem_vertex!, vertices, label_for, code_for
using MetaGraphsNext: inneighbors, outneighbors, has_edge, edge_labels, is_cyclic # TODO clean unused graph things
using MetaGraphsNext: DFSIterator
# 
# Overlays to composite
import ColorBlendModes
using ColorBlendModes: CompositeDestinationOver, BlendLighten, BlendMultiply
using Statistics: mean
# Vector graphics (this .svg is .xml)
import EzXML
using EzXML: readxml, write, Document, findfirst, setnodecontent!, firstnode
using EzXML: TextNode, AttributeNode, ElementNode, link!, unlink!
# User utilties
import Random
#
# Exports
#
export save_png_with_phys
export run_bitmapmap_pipeline, define_builder
# Export of pipeline functions
export establish_folder, unzip_tif, consolidate_elevation_data, water_overlay
export topo_relief, contour_lines_overlay, grid_overlay, ridge_overlay, summits_and_prominence
export join_layers, make_thumbnail_image, summits_regional_update, make_vector_graphics
# Export of builders and utilties
export SheetBuilder, SheetMatrixBuilder, full_folder_path
export northeast_corner, northeast_external_corner, northeast_internal_corner
export northwest_corner
export southeast_corner, southeast_external_corner, southeast_internal_corner
export southwest_corner, southwest_external_corner, southwest_internal_corner
export geo_grid_centre_single, geo_centre, bbox_external_string, polygon_string, show_augmented
export geo_area, nonzero_raster_closed_polygon_string, raster_polygon_string, nonzero_raster_rect
export copy_relevant_tifs_to_folder, tif_full_filenames_buried_in_folder, unzip_tif
export show_derived_properties, readclose
export utm32_to_33
# Export of manual edit tools:
export edit_in_gimp, open_as_temp_in_gimp, display_if_vscode

#
# Not exported because it's too generic: mark_at! and line!

const CONSOLIDATED_FNAM = "Consolidated.tif"
const TOPORELIEF_FNAM = "Toporelief.png"
const WATER_FNAM = "Water.png"
const CONTOUR_FNAM = "Contour.png"
const GRID_FNAM = "Grid.png"
const RIDGE_FNAM = "Ridge.png"
const MARKERS_FNAM = "Markers.png"
const THUMBNAIL_FNAM = "Thumbnail.png"
const SUMMITS_FNAM = "Summits.csv"
const COMPOSITE_FNAM = "Composite.png"
const PARSEABLE_FNAM = "Parse_builder.jl"
const ElevationGraph = MetaGraph{Int64, SimpleDiGraph{Int64}, Tuple{Int64, Int64}, Float32, Nothing, Set{@NamedTuple{min_x::Int64, min_y::Int64, max_x::Int64, max_y::Int64}}, typeof(identity), Nothing}

"""
TIFDIC :: Dict[String, NamedTuple}()

TIFDIC is a dictionary for storing 
non-zero boundary boxes for the .geo files. 

If files are being changed during runs, this is a 
potential source of errors. We consider it's worthwhile,
because  opening every file for every sheet is potentially
very time consuming and shouldn't be repeated unnecessarily.

On the other hand, the dictionary will be rebuilt and resaved if the file
is deleted.
"""
const TIFDIC = Dict{String, NamedTuple}()
const TIFDIC_FNAM = "tifdic.jls"

include("ini_file.jl")
include("png_phys.jl")
include("filter.jl")
include("palette.jl")
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
include("summits_on_sheet.jl")
include("summits_regional_poststep.jl")
include("vector_graphics.jl")
include("graph_utilties.jl")
include("graph_harvest.jl")
include("utilties_segmentation.jl")
include("utilties_edit_with_gimp.jl")


end