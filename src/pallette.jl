# Elevations map into three hypsometric colour pallettes.
# The pallettes are hard-coded so you must like them.

"""
    generate_directional_pallette_func()
    ---> generic function (elevation, direction_no::Int64) ---> RGB

These are the undiminished hue and luminosity. `relief_variable_exponent` modifies this colour based on terrain inclination.

The colours are based on 62°N early spring, a clear day with snow above a limit, time around 14:00.


`direction_no` == 1  - azimuth and elevation of the sun, nominally 202° at 9° above horizon)
`direction_no` == 2  - azimuth angle lateral to  the sun (nominally 322° at 30° above horizon). This is mainly the light from a clear blue sky.
`direction_no` == 3  - azimuth angle opposite of the sun (nominally 22° at 30° above horizon). This is the light thrown back from a clear blue sky.
`direction_no` == 4  - azimuth angle lateral to the sun (nominally 82° at 30° above horizon)

...where azimuth angle 0° is North and 90° is West.

# Example
```
julia> color_from_elevation_and_direction_no = generate_directional_pallette_func()
#11 (generic function with 1 method)

julia> color_from_elevation_and_direction_no(10, 2)
RGB{N0f8}(0.494,0.494,0.431)

```
"""
function generate_directional_pallette_func()
    elevation_limits = [
        1,   # Sea
        2,   # Foreshore
        3,   # Rocks by sea
        5,   # Early spring grass
        50,  # Early spring grass
        280, # Leaf forest without leaves
        400, # Part snow part rocks
        500] # Snow
    easter_map_sunny_colors = RGB{N0f8}.([
        RGB{Float64}(0.4633664202773951,0.6879114708933098,0.9262170302591296)
        RGB{Float64}(0.522,0.467,0.42)
        RGB{Float64}(0.722,0.663,0.608)
        RGB{Float64}(1.0,0.9325421705307276,0.6687703456742201)
        RGB{Float64}(0.888799950269018,0.8247253696065902,0.6444621877814288)
        RGB{Float64}(0.5179344009565391,0.5138638751188354,0.5331628857914726)
        RGB{Float64}(0.8211876549622903,0.8334831253595203,0.8530014255148916)
        RGB{Float64}(0.9880563911327263,0.9880556184197363,0.9961609305710946)])
    easter_map_shadow_colors = RGB{N0f8}.([
        RGB{Float64}(0.184,0.42,0.635)
        RGB{Float64}(0.173,0.192,0.255)
        RGB{Float64}(0.251,0.29,0.376)
        RGB{Float64}(0.32088349699205737,0.34628438217549173,0.34937541392838933)
        RGB{Float64}(0.4074036196362983,0.4244583520057037,0.3734053455910808)
        RGB{Float64}(0.25955638017789534,0.27526962361575374,0.33598587819752196)
        RGB{Float64}(0.30292319717824423,0.3673755209368712,0.5245856821866882)
        RGB{Float64}(0.47018902285535163,0.5784619137094618,0.6880473522836665)])
    # Side light is a mix, mostly shadow tints.
    # But snowy terrain looses detail from side light. By setting snowy terrain side colour
    # to black, snowy areas do not get any side light.
    easter_map_side_colors = 0.75 .* easter_map_shadow_colors .+ 0.25 * easter_map_sunny_colors
    easter_map_side_colors[8] = RGB{N0f8}(0.0,0.0,0.0)
    # Let's add some transitional colors between these:
    easter_500 = make_map_500(easter_map_sunny_colors, elevation_limits)
    easter_shadow_500 = make_map_500(easter_map_shadow_colors, elevation_limits)
    easter_side_500 = make_map_500(easter_map_side_colors,   elevation_limits)
    palette_matrix = hcat(easter_500, easter_side_500, easter_shadow_500, easter_side_500)
    return (elevation::Float32, direction_no::Int64) -> begin
        i = clamp(round(Int64, elevation), 1, 500)
        c = palette_matrix[i, direction_no]
        return c
    end
end


"""
    make_map_500(colors, upper_limits)

For hypsometric tints, locally adapted.
"""
function make_map_500(colors, upper_limits)
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
    RGB{N0f8}[color(z_above) for z_above in 1:500]
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
