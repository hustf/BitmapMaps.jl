# BitmapMaps.jl
Make topographic relief bitmaps for printing from 32-bit elevation data, overlain with vector graphics


# What does it do?

Make printable topographic relief maps. The foreground is vector graphics from [RouteMap.jl](https://github.com/hustf/RouteMap.jl),
and the background is topographic relief maps based on elevation data. The hypsometric colours resemble a clear, early afternoon in 
mid-February at 62°N, with snow cover above 500 m.

This example uses laser elevations with a 3 metre spacing. Publicly available data covers all of Norway with a spacing of 1 metre or less,
so we could potentially zoom in to the street level. Both surface (with trees and houses) and terrain (cleaned of such) are available.

<img src="resource/bitmap_detail.png" alt = "resource/bitmap_detail.png" style="display: inline-block; margin: 0 auto; max-width: 640px">

The example includes default elevation contours every 100 m with a fatter curve at 1000 m. The map projection,, unusually, follows the UTM grid, which is shown with 1 km spacing.

# How does it do it?

Rendering the finished bitmaps is not expected to run in one go. As a minimum, a folder hierarchy is established at the first go, user moves in data, and bitmaps are rendered after the second call.

Intermediate step results are saved in a folder structure defined in an .ini file. The .ini file is generated automatically, but can be updated by user.

If an intermediate file is deleted, the corresponding step is triggered to rerun.

Steps in the pipeline:

1) Define a `SheetMatrixBuilder` based on `home/BitmapMaps.ini`. Keywords overrule file values. Repl feedback for a preview of the bitmap's geographical extent and division into sheets (`define_builder`).
2) Establish folder hierarchy for storing intermediate data (`establish_folder`).
3) User action if missing: 
  - Make requests for public data at høydedata.no (requires email). Download to each sheet's folder.
  - Or, if data exists locally: (`copy_relevant_tifs_to_folder`)
4) Unzip elevation data (`unzip_tif`)
5) Consolidate elevation data (`consolidate_elevation_data`).
6) Identify water surfaces (`water_overlay`).
7) Make topographic reliefs (`topo_relief`) from a hardcoded hypsometric colour pallette.
8) Make elevation contours (`contour_lines_overlay`)
9) Make vector graphics and text covering the full map area. You may use RouteMap.jl for this step.
10) Make composite bitmaps: 
    - topographic reliefs 
    - water surfaces
    - elevation countours 
    - vector graphics and text

# Example
```
using BitmapMaps
smb = run_bitmapmap_pipeline(;complete_sheets_first) # This makes a default .ini and establishes all folders
copy_relevant_tifs_to_folder("yourpath", smb)
run_bitmapmap_pipeline()
```
# Current state

A first map has been made with scripting, ad hoc calculations and A4 sheets. We're streamlining the production of maps.

Code is being adapted from environment 'geoarrays', package 'RouteMap.jl' ' / example/ split ' and environment 'tutorial_images' / 'image segmentation.jl'.

We implemented and tested up to step 6 here, and have a skeleton for steps 7-8.

In the current commit, we struggled with finding good function names related to bounding boxes and (cropped extents). We still aren't satisfied. Does the extent
of a file include zero-padded regions? We have sort-of landed on using different, longer, function names for .tif files, loaded geoarrays, and builders.

Also currently, the type 'SheetMatrixBuilder' is over-determined. Using inches as a unit and then demanding integer values leads to some inconsistency. We'll get rid of that 
by changing the type - in the next revision?

Some of the changes from scripting workflow:

- Establish BitmapMaps.ini, tune default printer data (the scripted / manual workflow map was missing 1 mm due to 'random' variations during printing).
- Found a reliable way to print with actual scale. Use png's pHYs chunk, then print with an application that respects the settings. E.g. `Gimp`, `MS Paint`, and `IrFanview`.
- Identifying water surface is now done with a more advanced algorithm. Manual corrections will hopefully be less of a requirement.
- The colour palette is FixedPointNumbers.Normed{UInt8, 8}, not Float64 - based
- Introduce the SheetMatrixBuilder (iterator for printable sheets) and SheetBuilder (iterator for pixels in a sheet). Change sheet numbering to start in SW corner. See figure:

<img src="resource/matrix_sheet_cell_utm.svg" alt = "resource/matrix_sheet_cell_utm.svg" style="display: inline-block; margin: 0 auto; max-width: 640px">