# This file contains functions that are called from within the context of another sheet.

"""
    neighbour_folder(sb::SheetBuilder, direction::Symbol)
    ---> String

Returns an empty string if the neigbour folder does not exist.
"""
function neighbour_folder(sb::SheetBuilder, direction::Symbol)
    # Check argument
    dic_off = Dict( :n => (1, 0),
                    :s => (-1, 0),
                    :e => (0, 1),
                    :w => (0, -1))
    Δi, Δj = get(dic_off, direction, (0,0))
    Δi == Δj == 0 && throw(ArgumentError("Direction is $direction, ∉ [:n, :s, :e, :w]"))
    # Collect from sb
    fo = joinpath(splitpath(sb.pthsh)[1:end - 1])
    i = parse(Int, split(splitpath(sb.pthsh)[end])[1])
    j = parse(Int, split(splitpath(sb.pthsh)[end])[2])
    sw = southwest_external_corner(sb)
    ne = northeast_external_corner(sb)
    w = width_cell(sb) * cell_to_utm_factor(sb)
    h = height_cell(sb) * cell_to_utm_factor(sb)
    #
    fofo = ""
    fofo *= string(i + Δi) * " " * string(j + Δj)
    fofo *= "  " * string(sw[1] + Δj * w) * " " * string(sw[2] + Δi * h)
    fofo *= "  " * string(ne[1] + Δj * w) * " " * string(ne[2] + Δi * h)
    fullfo = joinpath(homedir(), fo, fofo)
    if isdir(fullfo)
        return fullfo
    else
        return ""
    end
end
function neighbour_folder_dict(sb)
    dic = Dict{Symbol, String}()
    for direction in [:s, :n, :w, :e]
        nefo = neighbour_folder(sb, direction)
        if nefo !== ""
            push!(dic, direction => nefo)
        end
    end
    dic
end

"""
    distribute_boundary_conditions(dic_neighbour, z, mea)

We create files in neigbouring folders, duplicating the edges of matrices
for elevation (z) and maximum elevation above (mea). 

The full matrices can be large and could be somewhat heavy to load. Therefore,
we store the edges (the boundary conditions for neighbouring sheets) as vectors
so as to save on memory. Still, calculating contact interaction requires loading two 
full matrices and up to eight  vectors.  
"""
function distribute_boundary_conditions(dic_neighbour, z, mea)
    dic_opposite = Dict(:n => :s, :e => :w, :s => :n, :w => :e)
    # If the dictionary is empty, nothing happens.
    if ! isempty(dic_neighbour)
        for (direction, destination_folder) in dic_neighbour
            # If we are, say, duplicating our northern boundary to
            # our neigbour folder in the north, that would be
            # the destination's south boundary.
            suffix = dic_opposite[direction]
            # Elevation
            fna_z = "boundary_z_$suffix.vec"
            save_boundary(joinpath(destination_folder, fna_z), z, direction)
            # Maximum elevation above
            fna_mea = "boundary_mea_$suffix.vec" 
            isdir(destination_folder) || throw(ArgumentError("The non-existing destination folder $(desitination_folder) should not have been put in dic_neighbour"))
            save_boundary(joinpath(destination_folder, fna_mea), mea, direction)
        end
    end
end
function save_boundary(ffna, matr, direction)
    if direction == :n
        out = matr[1, :]
    elseif direction == :s
        out = matr[end, :]
    elseif direction == :w
        out = matr[:, 1]
    elseif direction == :e
        out = matr[:, end]
    else
        throw(ArgumentError("direction"))
    end
    open(ffna, "w") do io
        writedlm(io, vec(out))
    end
end

function read_boundary_condition(fo, h_cell, w_cell)
    @assert isdir(fo) "$fo"
    z_s = read_bound(joinpath(fo, "boundary_z_s.vec"), w_cell)
    z_n = read_bound(joinpath(fo, "boundary_z_n.vec"), w_cell)
    z_w = read_bound(joinpath(fo, "boundary_z_w.vec"), h_cell)
    z_e = read_bound(joinpath(fo, "boundary_z_e.vec"), h_cell)
    mea_s = read_bound(joinpath(fo, "boundary_mea_s.vec"), w_cell)
    mea_n = read_bound(joinpath(fo, "boundary_mea_n.vec"), w_cell)
    mea_w = read_bound(joinpath(fo, "boundary_mea_w.vec"), h_cell)
    mea_e = read_bound(joinpath(fo, "boundary_mea_e.vec"), h_cell)
    z_s, z_n, z_w, z_e, mea_s, mea_n, mea_w, mea_e
end

function read_bound(ffna, n)
    if isfile(ffna)
        v = vec(readdlm(ffna, Float32))
        length(v) == n || throw("Unexpected boundary length. \n    $ffna \n    Expected: $n got: $(length(v))")
        typeof(v) <: Vector{Float32} || throw("Unexpected boundary type. \n    $ffna \n    Got: $(typeof(v)) \n    Size: $(size(v))")
    else
        v = zeros(Float32, n)
    end
    v
end

boundary_cond_zero(z) = boundary_cond_zero(eltype(z), size(z)...)
function boundary_cond_zero(numtype, h_cell, w_cell)
    horiz = zeros(numtype, w_cell)
    vert = zeros(numtype, h_cell)
    Tuple(repeat([horiz, horiz, vert, vert], 2))
end


function pick_mea_including_from_boundary(z_boundary, mea_boundary, z_i, elev_above)
    if z_boundary >= z_i
        # The boundary lies above, and may thus carry info about a taller summit outside of
        # our boundary.
        if mea_boundary > elev_above
            mea_boundary
        else
            elev_above
        end
    else
        # The boundary lies below this cell in the sheet, so if anything, the next sheet
        # is partly in this summit's region.
        elev_above
    end
end



function func_mea_contact!(maxtree, z, bcond)
    h_cell = size(z, 1)
    w_cell = size(z, 2)
    z_s, z_n, z_w, z_e, mea_s, mea_n, mea_w, mea_e = bcond
    to_cartesian(i) = Tuple(CartesianIndices((h_cell, w_cell))[i])
    @assert length(z_s) == w_cell
    @assert length(z_n) == w_cell
    @assert length(z_w) == h_cell
    @assert length(z_e) == h_cell
    @assert length(mea_s) == w_cell
    @assert length(mea_n) == w_cell
    @assert length(mea_w) == h_cell
    @assert length(mea_e) == h_cell

    # We're making a double closure below. This makes stack traces hard to read.
    function pick_relevant_value(i, elev_above)
        m, n = to_cartesian(i)
        if m == 1
            pick_mea_including_from_boundary(z_n[n], mea_n[n], z[i], elev_above)
        elseif m == h_cell
            pick_mea_including_from_boundary(z_s[n], mea_s[n], z[i], elev_above)
        elseif n == 1
            pick_mea_including_from_boundary(z_w[m], mea_w[m], z[i], elev_above)
        elseif n == w_cell
            pick_mea_including_from_boundary(z_e[m], mea_e[m], z[i], elev_above)
        else
            elev_above
        end
    end
    function recurse!(mea, elev_above, i)
        ea = pick_relevant_value(i, elev_above)
        previous_max_elevation_above = mea[i]
        if isnan(previous_max_elevation_above)
            # Remember this.
            mea[i] = ea
        elseif previous_max_elevation_above >= ea
            # Just leave quietly. Stop the recursion here.
            return mea
        else 
            # Forget about what we thought was the summit elevation.
            # Remember this instead!
            mea[i] = ea
        end
        # Let's look for our parent.
        parent_i = maxtree.parentindices[i]
        if parent_i == i
            # ⬆ We are our own parent, the root, lowest level in the sheet.
            # Stop the recursion here.
            return mea
        end
        # ⬇ Recurse: Tell current parent (which is lower down) 
        # about the leaf's elevation. It's the actual summit above us for all we know.
        recurse!(mea, ea, parent_i)
    end
end
