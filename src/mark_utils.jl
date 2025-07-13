# Lightweight utilities used for drawing lines and
# geometric symbols on existing bitmap images.
# Primarily used for drawing summit markers,
# in steps `summits_on_sheet` and `summits_regional_poststep`.
#

"""
    mark_at!(img, inds::AbstractArray{CartesianIndex{2}}, side, shapename::String)
    mark_at!(img, inds::AbstractArray{CartesianIndex{2}}; side::Int = 3, f_is_filled = func_is_on_square)
    mark_at!(img, I::CartesianIndex{2}; side::Int = 3, f_is_filled = func_is_on_square)

Mark in-place
- a shape centered on cartesian index 'I'.
- identical shapes centered on a collection 'inds'.

# Arguments (including keywords)

- side       defines the shape's square bounding box
- shapename  "on_square", "in_circle", or corresponding from `f_is_filled` below
- f_is_filled is a function generator taking the 'side' argument.
  The generated function take a tuple of indices and returns a Bool.

Predefined (not exported):
```
func_is_on_square
func_is_on_triangle
func_is_on_circle
func_is_in_square
func_is_in_triangle
func_is_in_circle
func_is_on_cross
func_is_on_xcross
func_is_on_hline
func_is_on_vline
```
"""
function mark_at!(img, inds::AbstractArray{CartesianIndex{2}}, side, shapename::String)
    dic = Dict(
       "on_square" => func_is_on_square,
       "on_triangle" => func_is_on_triangle,
       "on_circle" => func_is_on_circle,
       "in_square" => func_is_in_square,
       "in_triangle" => func_is_in_triangle,
       "in_circle" => func_is_in_circle,
       "on_cross" => func_is_on_cross,
       "on_xcross" => func_is_on_xcross,
       "on_hline" => func_is_on_hline,
       "on_vline" => func_is_on_vline)
    shapename ∈ keys(dic) || throw(ArgumentError("shapename $shapename"))
    f_is_filled = dic[shapename]
    mark_at!(img, inds; side, f_is_filled)
end

function mark_at!(img, inds::AbstractArray{CartesianIndex{2}}; side::Int = 3, f_is_filled = func_is_on_square)
    @assert isodd(side)
    maxoff = side ÷ 2
    Ω = displacement_offset(maxoff, f_is_filled)
    for I in inds
        _mark_at!(img, I, Ω)
    end
    img
end
function mark_at!(img, I::CartesianIndex{2}; side::Int = 3, f_is_filled = func_is_on_square)
    @assert isodd(side)
    maxoff = side ÷ 2
    Ω = displacement_offset(maxoff, f_is_filled)
    _mark_at!(img, I, Ω)
end

"""
    line!(img, A::CartesianIndex{2}, B::CartesianIndex{2}, tol_dist=sqrt(2)/2)

Draws a line from A to B in the given image, with a thickness defined by `tol_dist`.
"""
function line!(img, A::CartesianIndex{2}, B::CartesianIndex{2}, tol_dist=√2 / 2)
    i1, j1 = Tuple(A)
    i2, j2 = Tuple(B)
    li = abs(i2 - i1)
    lj = abs(j2 - j1)
    si = sign(i2 - i1)
    sj = sign(j2 - j1)
    max_dist2 = tol_dist^2
    # Get all valid indices for the image
    R = CartesianIndices(img)
    if li > lj
        err = li ÷ 2
        for _ in 1:(li + 1)
            color_neighbors!(img, R, CartesianIndex(i1, j1), max_dist2)
            err -= lj
            if err < 0
                j1 += sj
                err += li
            end
            i1 += si
        end
    else
        err = lj ÷ 2
        for _ in 1:(lj + 1)
            color_neighbors!(img, R, CartesianIndex(i1, j1), max_dist2)
            err -= li
            if err < 0
                i1 += si
                err += lj
            end
            j1 += sj
        end
    end
    img
end

"""
    color_neighbors!(img, R, C, max_dist2)

Colors all neighbors of index `C` within the squared distance `max_dist2` and bounds `R`.
"""
function color_neighbors!(img, R, C::CartesianIndex{2}, max_dist2)
    for di in -1:1
        for dj in -1:1
            neighbor = CartesianIndex(C[1] + di, C[2] + dj)
            if neighbor in R && (di^2 + dj^2 <= max_dist2)
                img[neighbor] = oneunit(eltype(img))
            end
        end
    end
    img
end


###########
# Internals
###########
function _mark_at!(img, p::CartesianIndex{2}, Ω)
    R = CartesianIndices(img)
    Ωₚ = filter(q -> in(q, R), Ref(p) .+ Ω)
    img[Ωₚ] .= oneunit(eltype(img))
    img
end

"""
    displacement_offset(max_offset, f_is_filled)
"""
function displacement_offset(max_offset, f_is_filled)
    rng = -max_offset:1:max_offset
    f = f_is_filled(max_offset)
    [CartesianIndex((i,j)) for i in rng, j in rng if f(i, j)]
end

function func_is_on_line(A::T, B::T) where T <: Tuple{Int64, Int64}
    s = hypot((A .- B)...)
    # This 'length' is the number of stops, not the distance:
    length = Int(ceil(s)) + 1
    ABi = range(A[1], B[1]; length)
    ABj = range(A[2], B[2]; length)
    tol_dist = √2 / 2
    max_dist2 = tol_dist^2
    (i, j) -> begin
        range_offsets_i = ABi .- i
        range_offsets_j = ABj .- j
        range_offsets = zip(range_offsets_i, range_offsets_j)
        any(range_offsets) do (di, dj)
            di^2 + dj^2 <= max_dist2
        end
    end
end

###################
# Shape definitions
##################

func_is_on_square(maxoff) = (i, j) -> abs(i) == maxoff || abs(j) == maxoff

function func_is_on_triangle(maxoff)
    # Circumradius from c.o.a.
    r = maxoff
    # End indices of base line
    aj = Int(round(r * √3 / 2))
    ai = Int(round(r /2))
    A = (-r, 0)
    B = (ai, aj)
    C = (ai, -aj)
    fAB = func_is_on_line(A, B)
    fBC = func_is_on_line(B, C)
    fCA = func_is_on_line(C, A)
    (i, j) -> fAB(i, j) || fBC(i, j) || fCA(i, j)
end

function func_is_on_circle(maxoff)
    radius = maxoff
    tol_dist = 0.5
    min_r2 = Int(round((radius - tol_dist)^2))
    max_r2 = Int(round((radius + tol_dist)^2))
    (i, j) -> begin
        min_r2 <= i^2 + j^2 <= max_r2
    end
end

func_is_on_hline(maxoff) =  (i, j) -> i == 0
func_is_on_vline(maxoff) =  (i, j) -> j == 0
function func_is_on_cross(maxoff)
    fh = func_is_on_hline(maxoff)
    fv = func_is_on_vline(maxoff)
    (i, j) -> fh(i, j) || fv(i, j)
end

function func_is_on_xcross(maxoff)
    r = maxoff
    A = (-r, -r)
    B = (r, r)
    C = (-r, r)
    D = (r, -r)
    fAB = func_is_on_line(A, B)
    fCD = func_is_on_line(C, D)
    (i, j) -> fAB(i, j) || fCD(i, j)
end

function func_is_in_triangle(maxoff)
    # Circumradius from c.o.a.
    r = maxoff
    # End indices of base line
    aj = Int(round(r * √3 / 2))
    ai = Int(round(r /2))
    A = (-r, 0)
    B = (ai, aj)
    C = (ai, -aj)
    function barycentric_coordinates(P)
        x1, y1 = A
        x2, y2 = B
        x3, y3 = C
        x, y = P
        detT = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
        α = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / detT
        β  = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / detT
        γ = 1 - α - β
        α, β, γ
    end
    (i, j) -> begin
        α, β, γ = barycentric_coordinates((i, j))
        0 <= α <= 1 && 0 <= β <= 1 && 0 <= γ <= 1
    end
end

function func_is_in_circle(maxoff)
    radius = maxoff
    tol_dist = 0.5
    max_r2 = Int(round((radius + tol_dist)^2))
    (i, j) -> begin
        i^2 + j^2 <= max_r2
    end
end

func_is_in_square(maxoff) = (i, j) -> true
