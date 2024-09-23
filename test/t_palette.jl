using BitmapMaps
using Test
using BitmapMaps: func_directional_pallette, RGB, N0f8, *, mark_at!
using BitmapMaps: easter_map_shadow_colors, easter_map_sunny_colors
using BitmapMaps: XYZ # Nice for rendering because multiplication doesn't change chromaticity
using ColorTypes: xyY # Nice for adjusting lightness without affecting chromaticity

##########################
# Compare, adjust palettes
##########################
# Old version
p0 = XYZ{Float32}.([
        RGB{Float64}(0.4633664202773951,0.6879114708933098,0.9262170302591296)
        RGB{Float64}(0.522,0.467,0.42)
        RGB{Float64}(0.722,0.663,0.608)
        RGB{Float64}(1.0,0.9325421705307276,0.6687703456742201)
        RGB{Float64}(0.888799950269018,0.8247253696065902,0.6444621877814288)
        RGB{Float64}(0.5179344009565391,0.5138638751188354,0.5331628857914726)
        RGB{Float64}(0.8211876549622903,0.8334831253595203,0.8530014255148916)
        RGB{Float64}(0.9880563911327263,0.9880556184197363,0.9961609305710946)
        RGB{Float64}(1.0,0.9883187896289786,0.9221551886496328)])
# Add yellow-red tint at high elevations
p1 = xyY{Float32}.([
    RGB{Float64}(0.4633664202773951,0.6879114708933098,0.9262170302591296)
    RGB{Float64}(0.522,0.467,0.42)
    RGB{Float64}(0.722,0.663,0.608)
    RGB{Float64}(1.0,0.9325421705307276,0.6687703456742201)
    RGB{Float64}(0.888799950269018,0.8247253696065902,0.6444621877814288)
    RGB{Float64}(0.5179344009565391,0.5138638751188354,0.5331628857914726)
    RGB{Float64}(0.8211876549622903,0.8334831253595203,0.8530014255148916)
    RGB{Float64}(0.9880563911327263,0.9880556184197363,0.9961609305710946)
    RGB{N0f8}(1.0,0.98,0.671)])
# 

# Adjust the lightness without changing the chromaticity
p2 = XYZ[
    xyY{Float32}(0.2318886f0,0.24823461f0,  0.405648 )
    xyY{Float32}(0.34537512f0,0.35115513f0, 0.4 )
    xyY{Float32}(0.33906987f0,0.34807217f0, 0.5 )
    xyY{Float32}(0.3679895f0,0.3966278f0,   0.8535153 )
    xyY{Float32}(0.35864177f0,0.38077757f0, 0.65106004 )
    xyY{Float32}(0.30870277f0,0.3204642f0,  0.4 )
    xyY{Float32}(0.3070763f0,0.32339594f0,  0.6620548 )
    xyY{Float32}(0.31152794f0,0.3270413f0,  0.97471976 )
    xyY{Float32}(0.3640501f0,0.40713286f0,  0.96879524 )]

# Lighter grass (4 and 5), forest (6)
p3 = XYZ[
    xyY{Float32}(0.2618911f0,0.2801038f0,1.5771875f0)
    xyY{Float32}(0.34537512f0,0.35115513f0, 0.9 )
    xyY{Float32}(0.33906987f0,0.34807217f0, 1.0 )
    xyY{Float32}(0.3679895f0,0.3966278f0,   0.97 )
    xyY{Float32}(0.35864177f0,0.38077757f0, 0.85 )
    xyY{Float32}(0.30870277f0,0.3204642f0,  0.6 )
    xyY{Float32}(0.3070763f0,0.32339594f0,  0.6620548 )
    xyY{Float32}(0.31152794f0,0.3270413f0,  0.97471976 )
    xyY{Float32}(0.3640501f0,0.40713286f0,  0.96879524 )]

vcat(p2', p3')

l = getfield.(xyY.(easter_map_shadow_colors()), :Y)

# Old version
s0 =     xyY{Float32}.([
    RGB{Float64}(0.184,0.42,0.635)
    RGB{Float64}(0.173,0.192,0.255)
    RGB{Float64}(0.251,0.29,0.376)
    RGB{Float64}(0.32088349699205737,0.34628438217549173,0.34937541392838933)
    RGB{Float64}(0.4074036196362983,0.4244583520057037,0.3734053455910808)
    RGB{Float64}(0.25955638017789534,0.27526962361575374,0.33598587819752196)
    RGB{Float64}(0.30292319717824423,0.3673755209368712,0.5245856821866882)
    RGB{Float64}(0.47018902285535163,0.5784619137094618,0.6880473522836665)
    RGB{Float64}(0.47018902285535163,0.5784619137094618,0.6880473522836665)])

# Relative ligtness compared to sun, old version
rl = getfield.(s0, :Y) ./ getfield.(p1, :Y)
rl .* l

# We keep the original shadow colors, but set the ligtness related to the new sunny colors p2
s1 =    XYZ{Float32}[
    xyY{Float32}(0.206187f0,0.21883054f0,  0.13680996 )
    xyY{Float32}(0.26566327f0,0.2672308f0, 0.06461003 )
    xyY{Float32}(0.26149404f0,0.26768488f0,0.083264075 )
    xyY{Float32}(0.30117953f0,0.3269658f0, 0.09553036 )
    xyY{Float32}(0.32338607f0,0.3573748f0, 0.14518279 )
    xyY{Float32}(0.27881297f0,0.28372413f0,0.10878235 )
    xyY{Float32}(0.24319595f0,0.24127704f0,0.112848 )
    xyY{Float32}(0.26290756f0,0.2841108f0, 0.2813371 )
    xyY{Float32}(0.26290756f0,0.2841108f0, 0.29431581 )]

# Lighter sea, grass (4 and 5), forest (6)
s2 =    XYZ{Float32}[
    xyY{Float32}(0.2618911f0,0.2801038f0,1.5771875f0)
    xyY{Float32}(0.26566327f0,0.2672308f0, 0.20 )
    xyY{Float32}(0.26149404f0,0.26768488f0,0.25 )
    xyY{Float32}(0.30117953f0,0.3269658f0, 0.19 )
    xyY{Float32}(0.32338607f0,0.3573748f0, 0.25 )
    xyY{Float32}(0.27881297f0,0.28372413f0,0.20 )
    xyY{Float32}(0.24319595f0,0.24127704f0,0.112848 )
    xyY{Float32}(0.26290756f0,0.2841108f0, 0.2813371 )
    xyY{Float32}(0.26290756f0,0.2841108f0, 0.29431581 )]


vcat(s1', s2')

# Lets show the four pallettes.
# Direction 1 (light from 202°, south-south-west) at top left
# Direction 2 (from 322°, north-west) at top-right
# Direction 3 (from 22°, oposite light from sun) at bottom-left
# Direction 4 (from 82°, light from east-north-east) at bottom-right
foo = func_directional_pallette()
w = 3000
pic = XYZ{Float32}[XYZ(1,1,1) for i in 1:w, j in 1:w]

for r in 0:1, c in 0:1
    i0 = r * w ÷ 2
    j0 = c * w ÷ 2
    dirno = 1 + r * 2 + c
    for i in 1:(w ÷ 2), j in 1:(w ÷ 2)
        ii = 1 + -i + i0 + (w ÷ 2)
        jj = j + j0
        I = CartesianIndex((ii, jj))
        pic[I] = (j / (w ÷ 2)) * foo(Float32(i), dirno)
    end
end
pic