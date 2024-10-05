# Utilty functions for smoothing, used in 'contour.jl' and 'ridges.jl'
# Finite Impulse Response Filter


"    highpass_coefficients(n)"
function highpass_coefficients(n)
    @assert isodd(n) && n > 1
    # Even number of coefficents on one side of center
    m = div(n,  2)
    # Cut-off frequency
    ωc = π / 16
    rng = -m:1:-1
    h = - sin.(ωc * rng) ./ (π * rng)
    # Center coefficent
    push!(h, 1 - ωc / π)
    # Mirror coefficients
    append!(h, reverse(h[1:m]))
end

"    lowpass_coefficients(n)"
function lowpass_coefficients(n)
    @assert isodd(n) && n > 1
    # Even number of coefficents on one side of center
    m = div(n,  2)
    # Cut-off frequency
    ωc = π / 16
    rng = -m:1:-1
    h = sin.(ωc * rng) ./ (π * rng)
    # Center coefficent
    push!(h, ωc / π)
    # Mirror coefficients
    append!(h, reverse(h[1:m]))
end


"Blackman window coefficients, n is an odd number"
function blackman_coefficients(n)
    @assert isodd(n) && n > 1
    m = div(n,  2)
    freqs = (-m:1:m) * π / m
    0.42 .+ 0.5 * cos.( freqs) .+ 0.08 * cos.(2 * freqs)
end

function fir_lp_coefficients(n)
    c = lowpass_coefficients(n) .* blackman_coefficients(n)
    # In order to have zero constant-signal attenuation regardless of length,
    # we scale the coeffients to sum one.
    # Also make this an OffsetArray suitable for kernel use
    centered(c ./ sum(c))
end
function fir_hp_coefficients(n)
    c = highpass_coefficients(n) .* blackman_coefficients(n)
    # Also make this an OffsetArray suitable for kernel use
    centered(c)
end

function conv(signal::Vector{T}, kernel::Vector{T}) where T
    len_signal = length(signal)
    len_kernel = length(kernel)
    @assert len_signal === len_kernel
    len_result = len_signal + len_kernel - 1
    result = zeros(T, len_result)
    for i in 1:len_signal
        for j in 1:len_kernel
            result[i + j - 1] += signal[i] * kernel[j]
        end
    end
    start_idx = div(len_kernel, 2) + 1
    end_idx = start_idx + len_signal - 1
    result[start_idx:end_idx]
end

 
"""
    smoothed_surface_fir(z; w = 49)

Apply a FIR low-pass filter with Blackman window. The intention is to avoid 
phase-distortion and to work with something familiar.

w is window length, an odd number.
"""
function smoothed_surface_fir(z; w = 49)
    # Coefficients, including a Blackman window.
    c = Float32.(fir_lp_coefficients(w))
    # Elevations smoothed by lp-filter
    imfilter(z, (c, transpose(c)), FIRTiled())
end


"""
    roughness_of_surface_fir(z; w = 69)

Apply a FIR high-pass filter with Blackman window. The intention is to avoid 
phase-distortion and to work with something familiar.

w is window length, an odd number.
"""
function roughness_of_surface_fir(z; w = 69)
    # Coefficients, including a Blackman window.
    c = Float32.(fir_hp_coefficients(w))
    # Elevations where the smooth variation is removed and
    # the bumps remain.
    imfilter(z, (c, transpose(c)), FIRTiled())
end


"""
    bumpiness(z; w = 69, cutoff⁻ = 1.5, z_max = 680.0f0, cutoff⁺ =  4.6477094f0)
    ---> Matrix{Gray{Float32}} (shape and size like z)

Returns positive values defined by a high-pass filter with length w, Blackman window.
Values defined as non-interesting by the other keywords are set to zero. 

Used by `is_forest`.
"""
function bumpiness(z; w = 69, cutoff⁻ = 1.5, z_max = 680.0f0, cutoff⁺ =  4.6477094f0)
    isodd(w) || throw(ArgumentError("w must be odd, not $w"))
    # The local bumps (+ and -)
    r = roughness_of_surface_fir(z; w)
    r1 = r .* (z .< z_max)
    # Take absolute value, and drop the large values, which are most likely artifacts
    # or houses. 
    map(r1) do ρ
        mag = abs(ρ)
        mag < cutoff⁻ ? Gray{Float32}(0.0f0) : mag > cutoff⁺ ? Gray{Float32}(0.0f0) : Gray{Float32}(mag / cutoff⁺)
    end
end


"""
    is_forest(z; forest_cells_min= 9062, w = 69, cutoff⁻ = 1.5, cutoff⁺ =  4.6477094f0, z_max = 680.0f0)
    is_forest(g::GeoArray)
    ---> Matrix{Gray{Bool}}  (shape and size like z)

Used by 'ridges' and 'contours'.

# Arguments

z                 Elevations. The other default parameters assumes a grid spacing of 1 utm meter.
forest_cells_min  Remove smaller forests
w                 Window length for high-pass filter
cutoff⁻           Drop filtered values below
cutoff⁺           Drop filtered values above (often steep ridges or artifacts)
z_max             Drop forest above this elevation
"""
function is_forest(z; forest_cells_min = 9062, w = 69, cutoff⁻ = 1.5, cutoff⁺ =  4.6477094f0, z_max = 680.0f0)
    # Based on roughness, values Gray(0.0f0 to 1.0f0⁻). 
    # Not too much and not too little. Not too high above the ocean.
    b = bumpiness(z; w, cutoff⁻, z_max, cutoff⁺)
    # Forests are interspersed with high and low values. Smooth that out.
    d = imfilter(b, Kernel.gaussian(2))
    # Change to black and white (forest)
    e = map(β -> Gray{N0f8}(β > 0.0), d)
    # Build a Segmented Image, each forest has its own label.
    # Each forest must be large enough to exclude typical rough cliffs (long and slender shapes wihtout much area).
    # One (or more) of the segments cover the non-foresty areas.
    g = felzenszwalb(e, 1, forest_cells_min)
    # Convert segmented image to a normal image matrix
    map(i-> Gray{Bool}(round(segment_mean(g, i))), labels_map(g))
end
is_forest(g::GeoArray) =  is_forest(transpose(g.A[:, :, 1]))
