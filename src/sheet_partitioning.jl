# Based on Luxor.Partition, adapted so as to ensure that an integer number of sheets
# cover the entire area.
# My printer tends to sometimes drop 1 mm of graphics when printer. 
# Add some pixels to the final output image to compensate for that. But that
# does not affect partitioning.

"""
    p = ShPartition(bmm_pixel_width, bmm_pixel_height, sheet_width, sheet_height)

A ShPartition is an iterator that, for each iteration, returns a tuple of:

  - the `x`/`y` point of the center of each tile in a set of sheets that divide up a
    rectangular space such as a page into rows and columns (relative to current 0/0)

  - the number of the tile

`bmm_pixel_width` and `bmm_pixel_height` are the dimensions of the area to be tiled,
`sheet_width`/`sheet_height` are the dimensions of the sheets.

Tiler and ShPartition are similar:

  - ShPartition lets you specify the width and height of a sheet

  - Tiler lets you specify how many rows and columns of sheets you want, and a margin

  ```

    sheets = ShPartition(1200, 1200, 30, 30)
    for (pos, n) in sheets
    
    # the point pos is the center of the tile
    
    end
  ```

You can access the calculated tile width and height like this:

```
    sheets = ShPartition(1200, 1200, 30, 30)
    for (pos, n) in sheets
        ellipse(pos.x, pos.y, sheets.sheet_width, sheets.sheet_height, :fill)
    end
```


It's sometimes useful to know which row and column you're currently on:

    sheets.currentrow
    sheets.currentcol

should have that information for you.

Unless the entire 'bitmapmap' (the entire map or canvas) and 'sheet' 
dimensions are exact multiples of each other, the bordering sheets will  
extend outside of the 'bitmapmap' with background content.

TODO: Hence, we need to define a background color.
TODO: Input:
    1)  utm bounding box, pixels sheet width and height, and an integer utm_to_pixel_factor.
        Output is number of sheets.
    2) utm south-west, pixels sheet width and height, number of sheets m and n, and an integer utm_to_pixel_factor.

    Point is also used for utm coordinates elsewhere, though that is float.

    numbering: South-west to south-east, ending at north-east.
               This is nice, almost necessary because the source elevation grid is measured this way.
               Hence, we don't lose data at the north border.

    Output from each iteration:
    (m, n) bb_utm bb_pix

    pixel
"""
mutable struct ShPartition
    bmm_pixel_width::Int
    bmm_pixel_height::Int
    sheet_width::Int
    sheet_height::Int
    nrows::Int
    ncols::Int
    currentrow::Int
    currentcol::Int
    function ShPartition(bmm_pixel_width, bmm_pixel_height, sheet_width, sheet_height)
        ncols = bmm_pixel_width / sheet_width |> ceil |> Integer
        nrows = bmm_pixel_height / sheet_height |> ceil |> Integer
        if ncols < 1 || nrows < 1
            throw(error("partition(): not enough space for sheets that size"))
        end
        currentrow = 1
        currentcol = 1
        new(bmm_pixel_width, bmm_pixel_height, sheet_width, sheet_height, nrows, ncols, currentrow, currentcol)
    end
end

# The first iteration, returns first and second state.
function Base.iterate(pt::ShPartition)
    x = -(pt.bmm_pixel_width / 2) + (pt.sheet_width / 2)
    y = (pt.bmm_pixel_height / 2) - (pt.sheet_height / 2)
    sheetnumber = 1
    x1 = x + pt.sheet_width
    y1 = y
    if (x1 + pt.sheet_width / 2) > (pt.bmm_pixel_width / 2)
        y1 -= pt.sheet_height
        x1 = -(pt.bmm_pixel_width / 2) + (pt.sheet_width / 2)
    end
    pt.currentrow, pt.currentcol = (div(sheetnumber - 1, pt.ncols) + 1, mod1(sheetnumber, pt.ncols))
    return ((Point(x, y), sheetnumber), (Point(x1, y1), sheetnumber + 1))
end

# Other than the first iteration, returns this and the next state.
function Base.iterate(pt::ShPartition, state)
    if state[2] > (pt.nrows * pt.ncols)
        return nothing
    end
    x = state[1].x
    y = state[1].y
    sheetnumber = state[2]
    x1 = x + pt.sheet_width
    y1 = y
    if (x1 + pt.sheet_width / 2) > (pt.bmm_pixel_width / 2)
        y1 -= pt.sheet_height
        x1 = -(pt.bmm_pixel_width / 2) + (pt.sheet_width / 2)
    end
    pt.currentrow, pt.currentcol = (div(sheetnumber - 1, pt.ncols) + 1, mod1(sheetnumber, pt.ncols))
    return ((Point(x, y), sheetnumber), (Point(x1, y1), sheetnumber + 1))
end


function Base.length(pt::ShPartition)
    pt.nrows * pt.ncols
end


function Base.getindex(pt::ShPartition, i::Int)
    1 <= i <= pt.ncols * pt.nrows || throw(BoundsError(pt, i))
    xcoord = -pt.bmm_pixel_width / 2 + mod1(i, pt.ncols) * pt.sheet_width - pt.sheet_width / 2
    ycoord = -pt.bmm_pixel_height / 2 + (div(i - 1, pt.ncols) * pt.sheet_height) + pt.sheet_height / 2
    return (Point(xcoord, ycoord), i)
end