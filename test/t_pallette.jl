# This was used to optimize the speed of hypsometric pallette.
using BitmapMaps
using Test
using BitmapMaps: generate_directional_pallette_func, RGB, N0f8
foo = generate_directional_pallette_func()
@test foo(100.1, 2) == RGB{N0f8}(0.482,0.482,0.427)
#using BenchmarkTools
# 15.054 ms (800013 allocations: 13.35 MiB)
# @btime [foo(x, dno) for x  in 0:0.1:10000, dno in 1:4]
