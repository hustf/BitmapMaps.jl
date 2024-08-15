# The algorithm and parameter values is optimized for in environments/tutorial_images/image_segmentation.jl
# Since these parameter values are well optimized, we use those hardcoded here.
#
# Work with the limited data in /resource.

using Test
using BitmapMaps
using BitmapMaps: mapwindow

###################################################
# Preparation (much the same as in `test_pipeline`)
###################################################

pth = "BitmapMaps\\test_topo"
# cleanup
if ispath(joinpath(homedir(), pth))
    sleep(1) # prevent resource busy error....
    rm(joinpath(homedir(), pth), recursive = true)
end

# Work with limited data from /resource
# Unzip the files to a temporary folder, where the folder name does not provide relevant info.
tmpdir_topo_relief = mktempdir()

let
    zfi = "../resource/eksport_800912_20240512.zip"
    zipfi = abspath(joinpath(@__DIR__, zfi))
    dest = abspath(joinpath(tmpdir_topo_relief, splitdir(zipfi)[2]))
    if ! isfile(dest)
        cp(zipfi, dest)
    end
    println(dest)
end


# A zip file now exists in tmpdir_topo_relief, as if downloaded by user.
# Extract and inspect.
unzip_tif(tmpdir_topo_relief)
fna = first(tif_full_filenames_buried_in_folder(tmpdir_topo_relief))
@test nonzero_raster_closed_polygon_string(fna)== "(18294 6937563, 18449 6937563, 18449 6937717, 18294 6937717, 18294 6937563)"
@test BitmapMaps.polygon_string(fna) ==  "POLYGON ((18294 6937562, 18449 6937562, 18449 6937717, 18294 6937717, 18294 6937562),\n                   (18294 6937563, 18449 6937563, 18449 6937717, 18294 6937717, 18294 6937563))"
##################
# More preparation
##################
smb = let
    mapscale = 1 // 1000
    density_pt_m⁻¹ = Int(round(1 / mapscale))
    side_sn_m = 6937717 - 6937563
    side_we_m = 18449 - 18294
    side_sn_paper_mm = Int(round(1000 * side_sn_m * mapscale))
    side_we_paper_mm = Int(round(1000 * side_we_m * mapscale))
    SheetMatrixBuilder(
        (18294, 6937563), # southwest_corner::Tuple{Int, Int},
        (1,1), # nrc::Tuple{Int, Int},
        1,  # cell_to_utm_factor::Int,
        side_we_paper_mm, # sheet_width_mm::Int,
        side_sn_paper_mm, # sheet_height_mm::Int,
        density_pt_m⁻¹,# density_pt_m⁻¹::Int,
        pth# pth::String
    )
end
@test size(smb[1].cell_iter) == (BitmapMaps.sheet_height_cell(smb), BitmapMaps.sheet_width_cell(smb))

# Establish test folder hierarchy
@test BitmapMaps.establish_folder.(smb) == [true]

# Now copy the relevant file to the relevant directory
fnas = copy_relevant_tifs_to_folder(tmpdir_topo_relief, smb)
# Consolidate
@test BitmapMaps.consolidate_elevation_data.(smb) == [true]
let
    for sb in smb
        fna = joinpath(full_folder_path(sb), BitmapMaps.CONSOLIDATED_FNAM)
        # Get elevation matrix
        za = let
            z = readclose(fna)
            transpose(z.A[:, :, 1])
        end
        @test sum(za) > 0
    end
end

################
# Render to file
################
let
    sb = smb[1]
    outfnam = joinpath(full_folder_path(sb), BitmapMaps.TOPORELIEF_FNAM)
    if isfile(outfnam)
        rm(outfnam)
    end
    BitmapMaps.topo_relief(sb)
    @test isfile(outfnam)
    if isfile(outfnam)
        rm(outfnam)
    end
end

###############################
# Render to display (in vscode)
###############################
let
    sb = smb[1]
    cell2utm = BitmapMaps.cell_to_utm_factor(sb)
    consfnam = joinpath(full_folder_path(sb), BitmapMaps.CONSOLIDATED_FNAM)
    @assert isfile(consfnam)
    g = readclose(consfnam)
    # We're transposing the source data here, because
    # it makes it easier to reason about north, south, east west.
    za = transpose(g.A[:, :, 1])
    img = BitmapMaps.__topo_relief(za, sb.cell_iter, cell2utm)
    display_if_vscode(img)
end



################################################################
# Verify directions for `_topo_relief` -> 'func_render'
################################################################
let
    # Prepare
    sb = smb[1]
    cell2utm = BitmapMaps.cell_to_utm_factor(sb)
    ny, nx = size(sb.cell_iter)
    source_indices = (1:cell2utm:(ny * cell2utm), 1:cell2utm:(nx * cell2utm))
    @assert nx == 155 && ny == 154 # Wider East-West than North-South
    consfnam = joinpath(full_folder_path(sb), BitmapMaps.CONSOLIDATED_FNAM)
    @assert isfile(consfnam)
    g = readclose(consfnam)
    # We're transposing the source data here, because
    # it makes it easier to reason about north, south, east west.
    za = transpose(g.A[:, :, 1])
    display_if_vscode(za)
    function f_verif(M::Matrix)
        @assert size(M) == (3, 3)
        _, w, _, n, z, s, _, e, _ = M
        deriv_east_west = (e - w) / 2
        deriv_south_north = (n - s) / 2
        mag = sqrt(1 + deriv_south_north^2 + deriv_east_west^2)
        n_ew = -deriv_east_west / mag
        n_sn = -deriv_south_north / mag
        n_up =  1 / mag
        # Unit vector in a right-handed Euclidean coordinate system
        n_ew, n_sn, n_up
    end
    #
    # Flat landscape
    #
    fill!(za, 1.0f0)
    ver = mapwindow(f_verif, za, (3, 3), indices = source_indices)
    # The unit normal vector points straight up
    @test ver[1] == (0, 0, 1)
    @test all(x -> x == ver[1], ver)
    #
    # Landscape sloping up 45°, highest in the north
    #
    z_northslope  = map(sb.cell_iter) do I
        easting, northing = sb.f_I_to_utm(I)
        northing - southwest_external_corner(sb)[2]
    end
    display_if_vscode(z_northslope)
    ver = mapwindow(f_verif, z_northslope, (3, 3), indices = source_indices)
    # The unit normal vector points south and up.
    # The normal at the edges depends on the default argument to mapwindow: border = "replicate"
    @test sum(abs.(ver[2] .- (0, -0.5√2, 0.5√2))) < 1e-10
    @test all(x -> x == ver[2], ver[2:(end - 2), :])
    #
    # Landscape sloping up 45°, highest in the east
    #
    z_eastslope  = map(sb.cell_iter) do I
        easting, northing = sb.f_I_to_utm(I)
        easting - southwest_external_corner(sb)[1]
    end
    display_if_vscode(z_eastslope)
    ver = mapwindow(f_verif, z_eastslope, (3, 3), indices = source_indices)
    # The unit normal vector points west and up
    @test sum(abs.(ver[2, 2] .- (-0.5√2, 0, 0.5√2))) < 1e-10
    # The normal at the edges depends on the default argument to mapwindow: border = "replicate"
    @test all(x -> x == ver[2, 2], ver[:, 2:(end - 1)])
end



#########
# Cleanup
#########
rm(tmpdir_topo_relief, recursive = true)
let 
    fo = full_folder_path(smb)
    if ispath(fo)
        sleep(1) # prevent resource busy error....
        rm(fo, recursive = true)
    end
end
if ispath(full_folder_path(smb))
    sleep(1) # prevent resource busy error....
    rm(full_folder_path(smb), recursive = true)
end


function f_verif1(M::Matrix)
    @assert size(M) == (3, 3)
    _, w, _, n, z, s, _, e, _ = M
    deriv_south_north = (n - s) / 2
    deriv_east_west = (w - e) / 2
    mag = sqrt(1 + deriv_south_north^2 + deriv_east_west^2)
    n_sn = -deriv_south_north / mag
    n_ew = -deriv_east_west / mag
    n_up =  1 / mag
    # Unit vector in a right-handed Euclidean coordinate system
    n_ew, n_sn, n_up
end 
