# This was used to optimize the speed of hypsometric palette.
using BitmapMaps
using Test
using BitmapMaps: func_directional_pallette, RGB, N0f8, *
foo = func_directional_pallette()
@test foo(100.1f0, 2) == RGB{N0f8}(0.506,0.565,0.443)
@test foo(500.1f0, 1) == RGB{N0f8}(0.988,0.988,0.996)
@test foo(500.0f0, 1) !== foo(1000.0f0, 1)
@test foo(500.0f0, 2) == foo(1000.0f0, 2) 
@test foo(500.0f0, 3) !== foo(1000.0f0, 3)
@test foo(500.0f0, 4) == foo(1000.0f0, 4)

# Lets show the four pallettes.
# Direction 1 (light from 202°, south-south-west) at top left
# Direction 2 (from 322°, north-west) at top-right
# Direction 3 (from 22°, oposite light from sun) at bottom-left
# Direction 4 (from 82°, light from east-north-east) at bottom-right

w = 3000
pic = RGB{N0f8}[RGB(1,1,1) for i in 1:w, j in 1:w]
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

#using BenchmarkTools
# 14.992 ms (800013 allocations: 13.35 MiB)
# 
#@btime [foo(x, dno) for x  in 0.0f0:0.1f0:10000, dno in 1:4]

