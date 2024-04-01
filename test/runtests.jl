using Test
olddir = pwd()
cd(mktempdir())
@testset "Read value from init file" begin
    include("t_ini_file.jl")
end
@testset "Set and inspect png physical print info chunk" begin
    include("t_png_phys.jl")
end
@testset "Test sheet partitioning" begin
    include("t_sheet_partitioning.jl")
end

cd(olddir)
