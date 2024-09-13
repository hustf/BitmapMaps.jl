# BitmapMaps.jl
Make topographic relief bitmaps for printing from 32-bit elevation data, overlain with vector graphics


# What does it do?

Make printable topographic relief maps. The foreground is vector graphics from [RouteMap.jl](https://github.com/hustf/RouteMap.jl),
and the background is topographic relief maps based on elevation data. The hypsometric colours resemble a clear, early afternoon in 
mid-February at 62°N, with snow cover above 500 m.

This example from an earlier version uses laser elevations with a 3 metre spacing. Publicly available data covers all of Norway with a spacing of 1 metre or less,
so we could potentially zoom in to the street level. Both surface (with trees and houses) and terrain (cleaned of such) are available.

<img src="resource/bitmap_detail.png" alt = "resource/bitmap_detail.png" style="display: inline-block; margin: 0 auto; max-width: 640px">

The example includes default elevation contours every 100 m with a fatter curve at 1000 m. The grid is UTM33N 1 km grid.

# How does it do it?

Rendering the finished bitmaps is not expected to run in one step. As a minimum:
- Step 1: `run_bitmapmap_pipeline()` => A folder hierarchy is established.

If needed, the pipeline establishes an .ini file with default argument values. User updates arguments in the file (recommended) or by passing keywords to the pipeline.

- Step 2: User moves relevant elevation data into user's project folder. It may cover a larger area. This may help: `copy_relevant_tifs_to_folder(source_folder, destination_folder)`.
- Step 3: `run_bitmapmap_pipeline()` => User inspects rendered images.

If an intermediate image file is deleted, the corresponding step is triggered to rerun.

If the keyword `complete_sheets_first` = true, which is the default, all steps are completed before the next sheet starts processing. Steps in the pipeline are:

1) Define a `SheetMatrixBuilder` based on `home/BitmapMaps.ini`. Keywords overrule file values. Repl feedback for a preview of the bitmap's geographical extent and division into sheets (`define_builder`).
2) Establish folder hierarchy for storing intermediate data (`establish_folder`).
3) User action if missing: 
  - Make requests for public data at høydedata.no (requires email). Download to each sheet's folder.
  - Or, if data exists locally: `copy_relevant_tifs_to_folder`
4) Unzip elevation data (`unzip_tif`).
5) Consolidate elevation data (`consolidate_elevation_data`).
6) Identify water surfaces (`water_overlay`).
7) Make topographic reliefs (`topo_relief`) from a hardcoded hypsometric colour palette.
8) Make elevation contours (`contour_lines_overlay`).
9) Add UTM grid (`grid_overlay`).
10) Mark dieders and ridges (`ridge_overlay`).
11) Find the prominence of summits and draw markers (`summit_markers`). Also writes an .csv file with matched geographical names.
12) Make a composite bitmap (`join_layers`).
13) Establish an .svg file and stylesheet, including bitmap and text elements (`make_vector_graphics`).
 
# Example
```
using BitmapMaps
smb = run_bitmapmap_pipeline(;complete_sheets_first) # This makes a default .ini and establishes all folders
copy_relevant_tifs_to_folder("yourpath", smb)
run_bitmapmap_pipeline()
```
# Current state

Steps 1-13 is fully working

Version 0.1.4 makes slight improvements to elevation contour lines. User can now modify the minimum length of "countour dashes", and the algorithm scales better when zoomed out.

: summit text and lake elevations is included. Style can be modified in the .ini file.

## Step 13, vector graphics:
Version 0.1.3: summit and lakes text has been refined to render well in all browsers and close to sheet edges. Style can be modified in the .ini file, text contents can be modified in .csv files.

## Step 11, summits:

Summit identification by name is first implemented in v0.0.40. Names are fetched via 'Stadnamn.jl' and added as a column to 'summits.csv'. A future
modfication would include overriding with user defined names in the .ini file, as a dictionary of "utm location" => "users own name".

Summit prominences are calculated taking boundary conditions between sheets into account. Since all boundary conditions are not available until the neighbouring sheets are calculated, this introduces iteration. Hence some summits will be appointed a lower prominence after re-running the pipeline.  

☛ Look for [ Info: No changes 
    to determine if more iteration is needed to find exact prominence values for summits.

Tall power lines are sometimes marked as obscure summits. Version v0.0.40 introduces a configurable filter, 'Markers / Mininum stress level' for filtering out unwanted summits. Stress is also reported in 'Summits.csv'. If some actual summits are missing, try to set the minimum stress level to 0 and iterate.

## General
`GeoArrays.jl` has breaking changes in  version 0.9 (we currently pin to 0.8.5). 

Some nice to know:

- Metadata for printing is made and stored in .png files. Settings are respected by e.g. `Gimp`, `MS Paint`, and `IrFanview`.
- Water surfaces often require manual touch-up. Try doing sea-level touch up in 'Consolidated.tif'. High elevations can also be 
  shown (reduced to 256 intensities in Gimp) by adding a temporary 'multiplication layer' in `Gimp`.
- The pipeline can reuse `SheetMatrixBuilder` and modify it by steps, using keywords.
- The UTM grid is the correct one for the local utm zone, even though data is from a country-wide zone. We might in the future add local UTM grid coordinates to
  summits as a tooltip (the <title> element).
- Sheet numbering starts in the SW corner. See figure:

<img src="resource/matrix_sheet_cell_utm.svg" alt = "resource/matrix_sheet_cell_utm.svg" style="display: inline-block; margin: 0 auto; max-width: 640px">


# Wishlist

- Side-by-side overview of all sheets
- Decrease the variation of contour lines with scaling.
- Add increasing warm tint to snow on higher elevations.

# Bounding box functions

Bounding boxes have meaning for:
   - GeoArrays (this type is defined by `GeoArrays.jl`)
   - file names referring GeoArrays
   - SheetMatrixBuilders (this package's main type)
   - SheetBuilders (this package's main type)

If you're inspecting your own job definitions, you may only need `show_augmented(smb)`.

`show_derived_properties` shows the interesting properties for file names and other types.
You may find `polygon_string` or `bbox_external_string` more useful for optimizing placement.

Why not just use `GeoArrays.bbox` and `GeoArrays.bbox_overlap`?
   - Two adjacent map sheets shares a boundary (x_max1 == x_min2), but do not overlap. In `GeoArrays.bbox_overlap`, two such boxes do overlap, because x_max1 refers a cell and not it's right edge.
   - In this package, UTM coordinates are integers (because that resolution is liberally licensed for all of Norway, and because we use folder names corresponding to external boundary boxes). GeoArrays.jl uses floating point numbers.
   - A sheet in a map book is naturally defined by its boundary. Such a boundary does not change with cell resolution or data density.
   - Downloaded elevation files may be zero-padded. We are mostly interested in the non-zero geographical region.
   - Rasters aren't simply matrices. Word definitions and conventions come from various professions.
   - When working with online map tools, we like to paste Well-Known-Text polygons.
