# Used by topo_relief, a step in pipeline.
# Elevations map into three hypsometric colour palettes.
# The palettes are hard-coded so you must like them.

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
    easter_map_side_colors = 0.75 .* easter_map_shadow_colors() .+ 0.25 * easter_map_sunny_colors()
    # We want more contrast in the snowy regions: Drop side light sources
    easter_map_side_colors[end - 1: end] .= RGB{N0f8}(0.0,0.0,0.0)
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
        RGB{N0f8}(c)
    end
    RGB{N0f8}[color(z_above) for z_above in 1:1500]
end

function easter_map_sunny_colors()
    #=
    RGB{N0f8}.([
        RGB{Float64}(0.4633664202773951,0.6879114708933098,0.9262170302591296)
        RGB{Float64}(0.522,0.467,0.42)
        RGB{Float64}(0.722,0.663,0.608)
        RGB{Float64}(1.0,0.9325421705307276,0.6687703456742201)
        RGB{Float64}(0.888799950269018,0.8247253696065902,0.6444621877814288)
        RGB{Float64}(0.5179344009565391,0.5138638751188354,0.5331628857914726)
        RGB{Float64}(0.8211876549622903,0.8334831253595203,0.8530014255148916)
        RGB{Float64}(0.9880563911327263,0.9880556184197363,0.9961609305710946)
        RGB{Float64}(1.0,0.9883187896289786,0.9221551886496328)])
    =#
    XYZ{Float32}[
    XYZ{Float32}(0.3793626f0,0.405648f0,0.85171324f0)
    XYZ{Float32}(0.18923497f0,0.4f0,0.16624333f0)
    XYZ{Float32}(0.39871234f0,0.5f0,0.36804605f0)
    XYZ{Float32}(0.79166275f0,0.8535153f0,0.5082442f0)
    XYZ{Float32}(0.6142646f0,0.65106004f0,0.44445688f0)
    XYZ{Float32}(0.22075173f0,0.4f0,0.26547974f0)
    XYZ{Float32}(0.62741256f0,0.6620548f0,0.75789607f0)
    XYZ{Float32}(0.9284165f0,0.97471976f0,1.0766958f0)
    XYZ{Float32}(0.9104392f0,0.96879524f0,0.5248448f0)]
end

function easter_map_shadow_colors()
    #=
    RGB{N0f8}.([
        RGB{Float64}(0.184,0.42,0.635)
        RGB{Float64}(0.173,0.192,0.255)
        RGB{Float64}(0.251,0.29,0.376)
        RGB{Float64}(0.32088349699205737,0.34628438217549173,0.34937541392838933)
        RGB{Float64}(0.4074036196362983,0.4244583520057037,0.3734053455910808)
        RGB{Float64}(0.25955638017789534,0.27526962361575374,0.33598587819752196)
        RGB{Float64}(0.30292319717824423,0.3673755209368712,0.5245856821866882)
        RGB{Float64}(0.47018902285535163,0.5784619137094618,0.6880473522836665)
        RGB{Float64}(0.47018902285535163,0.5784619137094618,0.6880473522836665)])
    =#
    XYZ{Float32}[
        XYZ{Float32}(0.12949122f0,0.13726963f0,0.36142537f0)
        XYZ{Float32}(0.030908918f0,0.13f0,0.054381445f0)
        XYZ{Float32}(0.06673847f0,0.13f0,0.12031099f0)
        XYZ{Float32}(0.08772189f0,0.13f0,0.10819712f0)
        XYZ{Float32}(0.13136747f0,0.1449441f0,0.12929884f0)
        XYZ{Float32}(0.061162192f0,0.14f0,0.09678767f0)
        XYZ{Float32}(0.11365008f0,0.15f0,0.24132648f0)
        XYZ{Float32}(0.2607117f0,0.28266907f0,0.44631496f0)
        XYZ{Float32}(0.2607117f0,0.4f0,0.44631496f0)
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
