using Test

olddir = pwd()
cd(mktempdir())
@testset "Read value from init file" begin
    include("t_ini_file.jl")
end
@testset "Set and inspect png physical print info chunk" begin
    include("t_png_phys.jl")
end
@testset "Hypsometric pallette" begin
    include("t_pallette.jl")
end
@testset "Test builders" begin
    include("t_builders.jl")
end
@testset "Test sheet builder utilties" begin
    include("t_builders_utilties.jl")
end
@testset "Test consolidate elevation data" begin
    include("t_unzip_and_consolidate.jl")
end
@testset "Test topo relief" begin
    include("t_topo_relief.jl")
end
@testset "Test identifying water surfaces" begin
    include("t_water.jl")
end
@testset "Test user utilties" begin
    include("t_user_geoarray_utilties.jl")
end
@testset "Test filter" begin
    include("t_filter.jl")
end
if ispath(joinpath(homedir(), "BitmapMaps\\Hov Litlevatn\\1 1  29000 6826535  31828 6830608"))
    @testset "Test elevation contour" begin # depends on local data files
        include("t_contour.jl")
    end
end
@testset "Test markers" begin
    include("t_markers.jl")
end
if ispath(joinpath(homedir(), "BitmapMaps", "proj 47675 6929520 57224 6947852", "1 1  47675 6929520  57224 6943269"))
    @testset "Test summit prominence" begin # depends on local data files
        include("t_summit_markers.jl")
    end
end
@testset "Test BitmapMaps pipeline" begin
    include("t_pipeline.jl")
end

cd(olddir)
