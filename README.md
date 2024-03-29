# BitmapMaps.jl
Make topographic relief bitmaps for printing from 32-bit elevation data, overlain with vector graphics


# What does it do?

Make printable topographic relief maps. The foreground is vector graphics from [RouteMap.jl](https://github.com/hustf/RouteMap.jl),
and the background is topographic relief maps based on elevation data. The hypsometric colours resemble a clear, early afternoon in 
mid-February at 62°N, with snow cover above 500 m.

This example uses laser elevations with a 3 metre spacing. Publicly available data covers all of Norway with a spacing of just 1 metre.
We could potentially zoom in to the street level.

<img src="resource/bitmap_detail.png" alt = "resource/bitmap_detail.png" style="display: inline-block; margin: 0 auto; max-width: 640px">

The example includes default elevation contours every 100 m with a fatter curve at 1000 m. The map projection is, unusually, follows the UTM grid, which is shown with 1 km spacing.

# How does it do it?

Rendering the finished bitmaps is not expected to run in one go. Intermediate step results are saved in a folder structure defined in an .ini file. The .ini file is generated automatically, but can be updated by user.

If an intermediate file is deleted, the corresponding step is triggered to rerun. 
Steps include:

- Make a grid for individual printed pages. Each page is associated with an utm coordinate bounding box.
- Establish folder hierarcy for storing intermediate data.
- Make requests for public data at høydedata.no (requires email). Download and unzip to appropriate folders.
- Sample and serialize elevation data.
- Identify water surfaces.
- Make topographic reliefs
- Make elevation contours
- Make vector graphics and text covering the full map area. You may use RouteMap.jl for this step.
- Make the composite bitmaps: 
    1) topographic reliefs 
    2) Water surfaces
    3) elevation countours 
    4) vector graphics and text


# Current state
Code currently resides in environment 'geoarrays' and in package 'RouteMap.jl' ' / example/ split '
