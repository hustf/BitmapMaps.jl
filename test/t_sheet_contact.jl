using BitmapMaps
using BitmapMaps: define_builder, neighbor_folder, neighbor_folder_dict, cell_to_utm_factor, full_folder_path
using BitmapMaps: read_boundary_condition, CONSOLIDATED_FNAM, MaxTree, distinct_summit_indices, maximum_elevation_above
using Test
smb = BitmapMaps.define_builder(;pth = "BitmapMaps/test")
sb = smb[1]
@test_throws ArgumentError neighbor_folder(sb, :nn)
@test neighbor_folder(sb, :n) == ""
@test neighbor_folder(sb, :s) == ""
@test neighbor_folder(sb, :e) == ""
@test neighbor_folder(sb, :w) == ""

fofo = joinpath(homedir(), "BitmapMaps", "proj 47675 6929520 57224 6947852", "1 1  47675 6929520  50858 6934103")
if ispath(fofo)

#= The following keywords should replicate the the following definitions in the .ini file:

[Printer consistent capability]
Printable width mm=191 # :sheet_width_mm
  #    Measured 193 mm. Allowing 2 mm  random variation.
Printable height mm=275 # :sheet_height_mm
  #    Measured 277 mm. Allowing 2 mm for random variation.
Stated density limit, dots per inch=600 # :density_limit_pt_inch⁻¹
  #    As advertised by Brother

[Geographical area]
Output paper sheets (rows columns)=(3 2) # :nrc
Southwest corner (utm easting northing)=(47675 6929520) # :southwest_corner
Output density, i.e. 'cells' / 'dots' / 'points' or 'pixels' per paper meter=16667 # :density_pt_m⁻¹
  #    For reference, 300  / inch = 300 / (0.0254 m) = 11811 m⁻¹ 
  #    Use lower values to slightly 'zoom in', otherwise use 'cell_to_utm_factor'.
Cell to utm factor, i.e. utm unit distance between elevation sampling points=2 # :cell_to_utm_factor
  #    How many 'utm metres' does a 'cell' / 'dot' / 'point' or 'pixel' side represent?
[File folder]
Top folders path under homedir()=BitmapMaps/proj 47675 6929520 57224 6947852 # :pth 
=#

# Establish 3x3 sheet folder structure (or make a full iteration)
smb = define_builder(cell_to_utm_factor = 1, nrc = (3, 3), southwest_corner = (47675, 6929520), 
    density_pt_m⁻¹ = 16667, pth = "BitmapMaps/proj 47675 6929520 57224 6947852")
@test smb.southwest_corner == (47675, 6929520)
smb = run_bitmapmap_pipeline(smb)


sb = smb[2, 2]
@test splitpath(neighbor_folder(sb, :n))[end] == "3 2  50858 6938686  54041 6943269"
@test splitpath(neighbor_folder(sb, :s))[end] == "1 2  50858 6929520  54041 6934103"
@test splitpath(neighbor_folder(sb, :e))[end] == "2 3  54041 6934103  57224 6938686"
@test splitpath(neighbor_folder(sb, :w))[end] == "2 1  47675 6934103  50858 6938686"
@test ! isempty(neighbor_folder_dict(sb))

# Let's have a look at mea...
sb = smb[3, 3]
ny, nx = size(sb.cell_iter)
si = CartesianIndices((1:cell_to_utm_factor(sb):(nx  * cell_to_utm_factor(sb)), 1:cell_to_utm_factor(sb):(ny * cell_to_utm_factor(sb))))
fofo = full_folder_path(sb)
@test isfile(joinpath(fofo, CONSOLIDATED_FNAM))
g = readclose(joinpath(fofo, CONSOLIDATED_FNAM))
z = transpose(g.A[si])
# Retrieve boundary conditions
bcond = read_boundary_condition(fofo, size(z, 1), size(z, 2))
maxtree = MaxTree(round.(z))
# Find the most important indices
summit_indices = distinct_summit_indices(z, maxtree)
  # Retrieve boundary conditions.
  # This is a tuple of four vectors (from file, or zero-filled with the correct length)
mea = maximum_elevation_above(z, bcond; maxtree, summit_indices)
imea = Int.(round.(mea))
if @isdefined levcols # levcols defined in `t_summits_on_sheet`
  indimg = levcols[imea]
end
# This value depends on how many times we have iterated. In the first run, we see the latter value.
@test abs(mea[1] - 1432.3724f0) < 0.01 || abs(mea[1] - 746.0909f0f0) < 0.01 

end # if ispath(fofo)
