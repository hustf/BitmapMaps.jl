using Test
using BitmapMaps
using BitmapMaps: ShPartition

shpa = ShPartition(3, 4, 1, 2)
@test length(shpa) == 6

(item, state) = iterate(shpa) # 1, 2
(item, state) = iterate(shpa, state) # 2,3
(item, state) = iterate(shpa, state) # 3,4
(item, state) = iterate(shpa, state) # 4, 5
(item, state) = iterate(shpa, state) # 5, 6
(item, state) = iterate(shpa, state) # 6, 7
next = iterate(shpa, state)
@test isnothing(next)

for (cenpos, n) in shpa
    @show n shpa.currentrow shpa.currentcol cenpos.x cenpos.y shpa.sheet_width shpa.sheet_height
    n == 1 && @test cenpos.x == -1
    n == 1 && @test cenpos.y == 1
    n == 6 && @test cenpos.x == 1
    n == 6 && @test cenpos.y == -1
    println( )
end


for (cenpos, n) in shpa
    @show n shpa.currentrow shpa.currentcol cenpos.x cenpos.y shpa.sheet_width shpa.sheet_height
    println( )
end

