# Step in pipeline.
# This renders a topographic relief map.

"""
    topo_relief(sb::SheetBuilder)
    topo_relief(fofo, cell_iter, cell2utm)
    ---> Bool

Output is an RGB colorspace .png image file.

Each rendered pixel is selected from four candidates, namely the lightest one.

Each candidate is calculated thus:

    Light source position & terrain normal vector => Lambertian reflection coefficient
    Elevation & lambertian reflection coefficient => reflection coefficient
    Light source colour * reflection coefficient => candidate colour & lightness

Note that the multiplication and lightness calculation is done in the XYZ colorspace.
The colours are defined in `pallette.jl`.
"""
function topo_relief(sb::SheetBuilder)
    topo_relief(full_folder_path(sb), sb.cell_iter, cell_to_utm_factor(sb))
end
function topo_relief(fofo, cell_iter, cell2utm)
    if isfile(joinpath(fofo, TOPORELIEF_FNAM))
        @debug "    $TOPORELIEF_FNAM in $fofo\n           already exists. Exiting `topo_relief`"
        return true
     end
    if ! isfile(joinpath(fofo, CONSOLIDATED_FNAM))
        @debug "    $CONSOLIDATED_FNAM in $fofo\n           does not exist. Exiting `topo_relief`"
        return false
    end
    res = RGB{N0f8}.(_topo_relief(fofo, cell_iter, cell2utm))
    # Feedback
    display_if_vscode(res)
    # Save
    ffna = joinpath(fofo, TOPORELIEF_FNAM)
    @debug "    Saving $ffna"
    save_png_with_phys(ffna, res)
    true
end
function _topo_relief(fofo, cell_iter, cell2utm)
    # Get elevation matrix. This samples every point regardless of cell_to_utm_factor
    za = elevation_full_res(fofo)
    __topo_relief(za, cell_iter, cell2utm)
end

function __topo_relief(za, cell_iter, cell2utm)
    # We need a 'render single output pixel function'. It changes a pixel at a time,
    # and its argument is its immediate surrounding in the source data.
    # It also needs to know more, but we're capturing that data.
    fr = func_render(func_directional_pallette())
    # Now map the render function to an output image.
    @debug "    Render topo relief"
    # Output image size
    ny, nx = size(cell_iter)
    # Source indices
    indices = (1:cell2utm:(ny * cell2utm), 1:cell2utm:(nx * cell2utm))
    # Apply
    mapwindow(fr, za, (3, 3); indices)
end

function func_render(f_hypso)
    # f_hypso takes two arguments: elevation and direction number 1..4.
    # It returns an XYZ{Float32}.
    (M::Matrix) -> @inbounds begin
        @assert size(M) == (3, 3)
        # If rows in M correspond to south -> north
        # and cols in M correspond to west -> east
        # _ n _
        # w z e
        # _ s _
        _, w, _, n, z, s, _, e, _ = M
        deriv_east_west = (e - w) / 2
        deriv_south_north = (n - s) / 2
        mag = sqrt(1 + deriv_south_north^2 + deriv_east_west^2)
        # Surface normal unit vectors, to the upper side.
        n_ew = -deriv_east_west / mag
        n_sn = -deriv_south_north / mag
        n_up =  1 / mag
        # Find the reflected color for light sources 1 to 4.
        vcol = reflected_color.(1:4, z, n_ew, n_sn, n_up, f_hypso)
        # We don't mix the light sources, but rather use the one
        # that is most luminous at this pixel.
        _, i = findmax(luminance, vcol)
        vcol[i]
    end
end

function reflection_coefficient(direction_no, z, n_ew, n_sn, n_up)
    @assert 1 <= direction_no <= 4
    sun = 202
    azimuth_deg = [sun, sun -180 -60, sun - 180, sun - 180 + 60]
    elevation_deg = [9, 30, 30, 30]
    az_deg = azimuth_deg[direction_no]
    el_deg = elevation_deg[direction_no]
    azim = az_deg * π / 180
    elev = el_deg * π / 180
    # Convert azimuth and elevation to a lighting direction vector.
    # Azimuth of 0     <=> Light from north, vector points north
    # Azimuth of π / 2 <=> Light from east, vector points east
    # Azimuth of π     <=> Light from south, vector points south
    # Elevation 0      <=> Light from horizon
    # Elevation π / 2  <=> Light from above
    #
    # Unit vector of the light source.
    l_ew = cos(elev) * sin(azim)
    l_sn = cos(elev) * cos(azim)
    l_up = sin(elev)
    # The dot product of light and surface normal is the fraction of light
    # reflected towards the observer
    lambert_reflection = lambert_shade(n_ew, n_sn, n_up, l_ew, l_sn, l_up)
    # We want a wider spread of reflected light further up, where snow is.
    # Note: lambert_reflection^shade_exponent(z) naively does not exceed 1.0. Thrust, but check.
    convert(Float32, min(1.0f0, lambert_reflection^shade_exponent(z, direction_no)))
end

function reflected_color(direction_no, z, n_ew, n_sn, n_up, f_hypso)
    r = reflection_coefficient(direction_no, z, n_ew, n_sn, n_up)
    r * f_hypso(z, direction_no)
end
function shade_exponent(z, direction_no)
    # Elevation dependent exponent for Lambertian shading.
    # We want the snow to return light from a wider cone.
    z1 = 400
    z2 = 500
    y1 = 1.2
    y2 = direction_no == 3 ? y1 : 0.4
    if z <= z1
        y1
    elseif z < z2
        linterp(y1, y2, z1, z2, z)
    else
        y2
    end
end
function lambert_shade(n_ew, n_sn, n_up, l_ew, l_sn, l_up)
    # Pure Lambertian reflection: dot product between surface normal and light direction normal.
    r = l_sn * n_sn + l_ew * n_ew + l_up * n_up
    # If the angle between l and n is more than 180°, reflection is negative. That means,
    # light would shine out from below the surface and could not be reflected from the surface.
    max(0.0, r)
end


function luminance(col::ColorTypes.RGB{N0f8})
    r, g, b = red(col), green(col), blue(col)
    # Linearization recommended by the recommended by the sRGB standard,
    r_lin = r <= 0.04045 ? r / 12.92 : ((r + 0.055) / 1.055) ^ 2.4
    g_lin = g <= 0.04045 ? g / 12.92 : ((g + 0.055) / 1.055) ^ 2.4
    b_lin = b <= 0.04045 ? b / 12.92 : ((b + 0.055) / 1.055) ^ 2.4
    # These coefficients reflect the human eye's sensitivity to different
    # wavelengths of light, as defined by the Rec. 709 standard.
    0.2126 * r_lin + 0.7152 * g_lin + 0.0722 * b_lin
end

luminance(col::ColorTypes.XYZ) = col.y

function lightness(col::ColorTypes.RGB{N0f8})
    lumi = luminance(col)
    # Calculate lightness from luminance
    lumi <= (216 / 24389) ? lumi * (24389 / 27) : lumi ^ (1 / 3) * 116 - 16
end
