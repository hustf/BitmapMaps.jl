# Used by topo_relief, a step in pipeline.
# Elevations map into three hypsometric colour palettes.
# The palettes and how they map to elevation are hard-coded.

"""
    func_directional_pallette()
    ---> generic function (elevation, direction_no::Int64) ---> RGB{N0f8}

These are the undiminished hue and luminosity. `relief_variable_exponent` modifies this colour based on terrain inclination.

The colours are based on 62°N early spring, a clear day with snow above a limit, time around 14:00.


`direction_no` == 1  - azimuth and elevation of the sun, nominally 202° at 9° above horizon)
`direction_no` == 2  - azimuth angle lateral to  the sun (nominally 322° at 30° above horizon). This is mainly the light from a clear blue sky.
`direction_no` == 3  - azimuth angle opposite of the sun (nominally 22° at 30° above horizon). This is the light thrown back from a clear blue sky.
`direction_no` == 4  - azimuth angle lateral to the sun (nominally 82° at 30° above horizon)

...where azimuth angle 0° is North and 90° is West.

# Example
```
julia> color_from_elevation_and_direction_no = func_directional_pallette()
#11 (generic function with 1 method)

julia> color_from_elevation_and_direction_no(10, 2)
RGB{N0f8}(0.494,0.494,0.431)

```
"""
function func_directional_pallette()
    elevation_limits = [
        1,   # Sea
        2,   # Foreshore
        3,   # Rocks by sea
        5,   # Early spring grass
        50,  # Early spring grass
        280, # Leaf forest without leaves
        400, # Part snow part rocks
        500, # Snow
       1500] # Snow, reddish tinted for elevation vis.
    # Side light is a mix, mostly shadow tints.
    # But snowy terrain looses detail from side light. By setting snowy terrain side colour
    # to black, snowy areas do not get any side light.
    easter_map_side_colors = 0.75f0 .* easter_map_shadow_colors() .+ 0.25f0 .* easter_map_sunny_colors()
    # We want more contrast in the snowy regions: Drop side light sources
    easter_map_side_colors[end - 1: end] .= XYZ{Float32}(0.0,0.0,0.0)
    # Let's add some transitional colors between these:
    easter_1500 = make_map_1500(easter_map_sunny_colors(), elevation_limits)
    easter_shadow_1500 = make_map_1500(easter_map_shadow_colors(), elevation_limits)
    easter_side_1500 = make_map_1500(easter_map_side_colors,   elevation_limits)
    palette_matrix = hcat(easter_1500, easter_side_1500, easter_shadow_1500, easter_side_1500)
    f = let palette_matrix = palette_matrix, max_elevation = maximum(elevation_limits)
        (elevation::Float32, direction_no::Int64) -> begin
            i = clamp(round(Int64, elevation), 1, max_elevation)
            palette_matrix[i, direction_no]
        end
    end
    f
end


"""
    make_map_1500(colors, upper_limits)

For hypsometric tints, locally adapted.
"""
function make_map_1500(colors, upper_limits)
    @assert length(colors) == length(upper_limits)
    function color(z_above)
        i = findlast(<=(z_above), upper_limits)
        if i < length(upper_limits)
            z1 = upper_limits[i]
            z2 = upper_limits[i + 1]
            c1 = colors[i]
            c2 = colors[i + 1]
            c = linterp(c1, c2, z1, z2, z_above)
        else
            c = colors[end]
        end
        c
    end
    XYZ{Float32}[color(z_above) for z_above in 1:1500]
end

function easter_map_sunny_colors()
    XYZ{Float32}[
        XYZ{Float32}(0.3226512f0,0.3441505f0,0.565619f0) # (1.4746368f0,1.5771875f0,2.5789008f0)
        XYZ{Float32}(0.88518596f0,0.9f0,0.77778375f0)
        XYZ{Float32}(0.9741367f0,1.0f0,0.8988308f0)
        XYZ{Float32}(0.8999617f0,0.97f0,0.575656f0)
        XYZ{Float32}(0.8005868f0,0.85f0,0.58168757f0)
        XYZ{Float32}(0.5779793f0,0.6f0,0.6943048f0)
        XYZ{Float32}(0.6286453f0,0.6620548f0,0.75649554f0)
        XYZ{Float32}(0.9284835f0,0.97471976f0,1.0772146f0)
        XYZ{Float32}(0.86627734f0,0.96879524f0,0.5444829f0)
        ]
end

function easter_map_shadow_colors()
    XYZ{Float32}[
        XYZ{Float32}(0.3226512f0,0.3441505f0,0.565619f0) # (1.4746368f0,1.5771875f0,2.5789008f0)
        XYZ{Float32}(0.19882682f0,0.2f0,0.34958988f0)
        XYZ{Float32}(0.24421817f0,0.25f0,0.4397158f0)
        XYZ{Float32}(0.17501558f0,0.19f0,0.21608494f0)
        XYZ{Float32}(0.22622333f0,0.25f0,0.22332238f0)
        XYZ{Float32}(0.19653808f0,0.2f0,0.308372f0)
        XYZ{Float32}(0.113745496f0,0.112848f0,0.2411178f0)
        XYZ{Float32}(0.26034087f0,0.2813371f0,0.44855928f0)
        XYZ{Float32}(0.27235097f0,0.29431581f0,0.46925232f0)
    ]
end


"""
    linterp(y1, y2, x1, x2, x3)
    linterp(y1, y2, normx)
    ---> y3 or normy

Linear interpolation.
"""
function linterp(y1, y2, x1, x2, x3)
    @assert x2 > x1
    normx = (x3 - x1) / (x2 - x1)
    linterp(y1, y2, normx)
end
function linterp(y1, y2, normx)
    @assert 0 <= normx <= 1
    y1 * (1 - normx) + y2 * normx
end
