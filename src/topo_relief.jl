# Step in pipeline.
# This renders a topographic relief map.

"""
    topo_relief(sb::SheetBuilder)
    topo_relief(fofo, cell_iter, cell2utm, f_hypso, f_reflect)
    ---> Bool

Output is an RGB colorspace .png image file.

Each rendered pixel is selected from four candidates, namely the lightest one.

Each candidate is calculated thus:

    Light source position & terrain normal vector => Lambertian reflection coefficient
    Elevation & lambertian reflection coefficient => reflection coefficient
    Light source colour * reflection coefficient => candidate colour & lightness

Note that the multiplication and lightness calculation is done in the XYZ colorspace.
The colours are defined in `pallette.jl`.

For customization of pallette, replace `f_hypso`. For customization of lighting, replace 
`f_reflect`. The argument list must be equivalent to `func_directional_pallette` and 
`func_reflection_coefficient`.
"""
function topo_relief(sb::SheetBuilder)
    # This function determines the color of the light source for a cell (a pixel in a raster)
    # The color varies with elevation above sea and light source number 
    f_hypso = func_directional_pallette()
    # This function determines how much of that color is reflected, given light source number,
    # the surface normal vector and elevation above sea.
    f_reflect = func_reflection_coefficient()
    #
    topo_relief(full_folder_path(sb), sb.cell_iter, cell_to_utm_factor(sb), f_hypso, f_reflect)
end
function topo_relief(fofo, cell_iter, cell2utm, f_hypso, f_reflect)
    if isfile(joinpath(fofo, TOPORELIEF_FNAM))
        @debug "    $TOPORELIEF_FNAM in $fofo\n           already exists. Exiting `topo_relief`"
        return true
     end
    if ! isfile(joinpath(fofo, CONSOLIDATED_FNAM))
        @debug "    $CONSOLIDATED_FNAM in $fofo\n           does not exist. Exiting `topo_relief`"
        return false
    end
    res = RGB{N0f8}.(_topo_relief(fofo, cell_iter, cell2utm, f_hypso, f_reflect))
    # Feedback
    display_if_vscode(res)
    # Save
    ffna = joinpath(fofo, TOPORELIEF_FNAM)
    @debug "    Saving $ffna"
    save_png_with_phys(ffna, res)
    true
end
function _topo_relief(fofo, cell_iter, cell2utm, f_hypso, f_reflect)
    # Get elevation matrix. This samples every point regardless of cell_to_utm_factor
    za = elevation_full_res(fofo)
    __topo_relief(za, cell_iter, cell2utm, f_hypso, f_reflect)
end

function __topo_relief(za, cell_iter, cell2utm, f_hypso, f_reflect)
    # We need a 'render single output pixel function'. It changes a pixel at a time,
    # and its argument is its immediate surrounding in the source data.
    # It also needs to know more, but we're capturing that data.
    fr = func_render(f_hypso, f_reflect)
    # Now map the render function to an output image.
    @debug "    Render topo relief"
    # Output image size
    ny, nx = size(cell_iter)
    # Source indices
    indices = (1:cell2utm:(ny * cell2utm), 1:cell2utm:(nx * cell2utm))
    # Apply
    mapwindow(fr, za, (3, 3); indices)
end

function func_render(f_hypso, f_reflect)
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
        deriv_west_east = (e - w) / 2
        deriv_south_north = (n - s) / 2
        mag = sqrt(1 + deriv_south_north^2 + deriv_west_east^2)
        # Surface normal unit vectors, to the upper side.
        n_we = -deriv_west_east / mag
        n_sn = -deriv_south_north / mag
        n_up =  1 / mag
        # Find the reflected color for light sources 1 to 4.
        vcol = reflected_color.(1:4, z, n_we, n_sn, n_up, f_hypso, f_reflect)
        # We don't mix the light sources, but rather use the one
        # that is most luminous at this pixel.
        _, i = findmax(luminance, vcol)
        vcol[i]
    end
end

"""
    reflected_color(direction_no, z, n_we, n_sn, n_up, f_hypso, f_reflect)
    ---> XYZ{Float32}
```
"""
function reflected_color(direction_no, z, n_we, n_sn, n_up, f_hypso, f_reflect)
    f_hypso(z, direction_no) * f_reflect(direction_no, z, n_we, n_sn, n_up)
end

"""
    func_reflection_coefficient(; sun_deg = 202, light_elev_deg = [9, 30, 30, 30])
    ---> generic function (direction_no, z, n_we, n_sn, n_up) ---> Float32
"""
function func_reflection_coefficient(; sun_deg = 202, light_elev_deg = [9, 30, 30, 30])
    # Direction no: From sun, opposite sun, one side 60°, other side 60 °
    Δazim_deg = [0, -180 - 60, - 180, - 180 + 60]
    light_azim_deg = sun_deg .+ Δazim_deg
    light_azim = light_azim_deg .* (π / 180)
    light_elev = light_elev_deg .* (π / 180)
    #
    # Convert light azimuth and light elevation to x-y-z light unit vectors.
    # Azimuth of 0     <=> Light from north, vector points north
    # Azimuth of π / 2 <=> Light from east, vector points east
    # Azimuth of π     <=> Light from south, vector points south
    # Elevation 0      <=> Light from horizon
    # Elevation π / 2  <=> Light from above
    #
    # Unit vector components of all the light sources
    l_ew = cos.(light_elev) .* sin.(light_azim)
    l_sn = cos.(light_elev) .* cos.(light_azim)
    l_up = sin.(light_elev)
    # Closure on the light source vectors
    f = let l_ew = l_ew, l_sn = l_sn, l_up = l_up
        (dno, z, n_we, n_sn, n_up) -> begin
            @assert 1 <= dno <= 4
            # The dot product of light and surface normal is the fraction of light
            # reflected towards the observer. This is 'lambert shade'.
            #
            # We modify the lambert shading where there is snow.
            # Lambert_reflection^shade_exponent(z, dno) naively does not exceed 1.0. 
            # Thrust, but check.
            min(1.0f0, convert(Float32, lambert_shade(n_we, n_sn, n_up, l_ew[dno], l_sn[dno], l_up[dno]) ^
                shade_exponent(z, dno)))
        end
    end
end


"""
    shade_exponent(z, direction_no)
    ---> Float64

Hard-coded values are:

Elevation     z <= 400 --> 1.2 
        400 < z < 500  --> interpolated
        500 <= z       --> 0.4 for light sources 1 to 3
                          1.2 for light source 3 (opposite sun)

When used as an exponent to Lambertian shading, this means:

Elevation     z <= 400 --> Slightly blank surface, reflects a little less all light angles considered
        400 < z < 500  --> Transition
        500 <= z       --> For light source 1-3: Matte appearance, reflects more  all light angles considered
                           For light source 3 (opposite sun), same for all elevations
"""
function shade_exponent(z, direction_no)
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

"""
    lambert_shade(n_we, n_sn, n_up, l_ew, l_sn, l_up)
    ---> Float64
"""
function lambert_shade(n_we, n_sn, n_up, l_ew, l_sn, l_up)
    # Pure Lambertian reflection: dot product between surface normal and light direction normal.
    r = l_sn * n_sn + l_ew * n_we + l_up * n_up
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
