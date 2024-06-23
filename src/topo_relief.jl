# This renders a topographic relief map.

"""
    topo_relief(sb::SheetBuilder)
    topo_relief(fofo)

"""
function topo_relief(sb::SheetBuilder)
    topo_relief(full_folder_path(sb), sb.cell_iter, cell_to_utm_factor(sb), sb.density_pt_m⁻¹)
end
function topo_relief(fofo, cell_iter, cell2utm, density_pt_m⁻¹)
    if isfile(joinpath(fofo, TOPORELIEF_FNAM))
        @debug "$TOPORELIEF_FNAM in $fofo already exists. Exiting `topo_relief`."
        return true
     end
    if ! isfile(joinpath(fofo, CONSOLIDATED_FNAM))
        @debug "$CONSOLIDATED_FNAM in $fofo does not exist. Exiting `topo_relief`."
        return false
    end
    res = _topo_relief(fofo, cell_iter, cell2utm, density_pt_m⁻¹)
    ffna = joinpath(fofo, TOPORELIEF_FNAM)
    @debug "Saving $ffna"
    save_png_with_phys(ffna, res, density_pt_m⁻¹)
    true
end
function _topo_relief(fofo, cell_iter, cell2utm, density_pt_m⁻¹)
    fr = generate_render_func(generate_directional_pallette_func())
    # Get elevation matrix. This samples every point regardless of cell_to_utm_factor
    g = readclose(joinpath(fofo, CONSOLIDATED_FNAM))
    za = g.A[:, :, 1]
    # Establish output image matrix
    ny, nx = size(cell_iter)
    source_indices = (1:cell2utm:(nx * cell2utm), 1:cell2utm:(ny * cell2utm))
    @debug "Render topo relief"
    relief = mapwindow(fr, za, (3, 3), indices = source_indices)
    display(transpose(relief))
    transpose(relief)
end

function generate_render_func(f_hypso)
    # f_hypso takes two arguments: elevation and direction number 1..4.
    # It returns an RGB{N0f8}.
    (M::Matrix) -> @inbounds begin
        @assert size(M) == (3, 3)
        _, w, _, n, z, s, _, e, _ = M
        deriv_south_north = (n - s) / 2
        deriv_east_west = (w - e) / 2
        # Surface normal unit vectors, to the upper side.
        mag = sqrt(1 + deriv_south_north^2 + deriv_east_west^2)
        n_sn = -deriv_south_north / mag
        n_ew = -deriv_east_west / mag
        n_up =  1 / mag
        # Find the color reflected from the sun's direction
        col  = reflected_color(1, z, n_sn, n_ew, n_up, f_hypso)
        for lightsource_no in 2:4
            othercol  = reflected_color(lightsource_no, z, n_sn, n_ew, n_up, f_hypso)
            if luminance(othercol) > luminance(col)
                col = othercol
            end
        end
        col
    end
end

function reflected_color(direction_no, z, n_sn, n_ew, n_up, f_hypso)
    @assert 1 <= direction_no <= 4
    expon = shade_exponent(z)
    sun = 202
    azimuth_deg = [sun, sun -180 -60, sun - 180, sun - 180 + 60]
    elevation_deg = [9, 30, 30, 30]
    az_deg = azimuth_deg[direction_no]
    el_deg = elevation_deg[direction_no]
    lambert_reflection = lambert_shade(n_sn, n_ew, n_up, az_deg * π / 180, el_deg * π / 180)
    reflection = convert(N0f8, min(1.0f0, lambert_reflection^expon))
    f_hypso(z, direction_no) * reflection
end
function shade_exponent(z)
    # Elevation dependent exponent for Lambertian shading.
    # We want the snow to return light from a wider cone.
    upper_limit_exponent_1 = 400
    lower_limit_low_exponent = 500
    low_exponent = 0.3
    if z <= upper_limit_exponent_1
        1.0
    elseif z >= lower_limit_low_exponent
        low_exponent
    else
        1 + (low_exponent - 1) * (z - upper_limit_exponent_1) / (lower_limit_low_exponent - upper_limit_exponent_1)
    end
end
function lambert_shade(n_sn, n_ew, n_up, azim, elev)
    # Convert azimuth and elevation to a lighting direction vector.
    # Azimuth of 0     <=> Light from north, vector points north
    # Azimuth of π / 2 <=> Light from east, vector points east
    # Azimuth of π     <=> Light from south, vector points south
    # Elevation 0      <=> Light from horizon
    # Elevation π / 2  <=> Light from above
    le = cos(elev) * sin(azim)
    ln = cos(elev) * cos(azim)
    lu = sin(elev)
    # Vector from the point in question to the light source.
    # Same coordinate system as surface normals (just to be clear about positive directions)
    l_sn = ln
    l_ew = le
    l_up = lu
    #
    # Pure Lambertian reflection: dot product between surface normal and light direction normal.
    r = l_sn * n_sn + l_ew * n_ew + l_up * n_up
    # If the angle between l and n is more than 180°, reflection is negative. That means,
    # light would shine out from below the surface and could not be reflected from the surface.
    max(0.0, r)
end


function luminance(col::ColorTypes.RGB{N0f8})
    r = col.r / 255
    g = col.g / 255
    b = col.b / 255
    r_lin = r <= 0.04045 ? r / 12.92 : ((r + 0.055) / 1.055) ^ 2.4
    g_lin = g <= 0.04045 ? g / 12.92 : ((g + 0.055) / 1.055) ^ 2.4
    b_lin = b <= 0.04045 ? b / 12.92 : ((b + 0.055) / 1.055) ^ 2.4
    0.2126 * r_lin + 0.7152 * g_lin + 0.0722 * b_lin
end

# For inspection during development (dead)
function to_gray(z)
   mi, ma = extrema(z)
   zr = (z .- mi) ./ (ma - mi)
   Gray{N0f8}.(zr)
end

function lightness(col::ColorTypes.RGB{N0f8})
    lumi = luminance(col)
    # Calculate lightness from luminance
    lumi <= (216 / 24389) ? lumi * (24389 / 27) : lumi ^ (1 / 3) * 116 - 16
end
