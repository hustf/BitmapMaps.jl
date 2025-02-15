# Utilty functions for smoothing, used in 'contour.jl' and 'ridges.jl'
# - Windowed Finite Impulse Response Filter
# - General terrain smoothing
# - Quantify roughness 
# - Group patches based on roughness

"    highpass_coefficients(n; nyquist_denom::Int = 2)"
function highpass_coefficients(n; nyquist_denom::Int = 2)
    isodd(n) && n > 1 || throw(ArgumentError("n"))
    nyquist_denom > 1 || throw(ArgumentError("nyquist_denom"))
    # Even number of coefficents on one side of center
    m = div(n,  2)
    # Cut-off frequency
    ωc = π / nyquist_denom
    rng = -m:1:-1
    h = - sin.(ωc * rng) ./ (π * rng)
    # Center coefficent
    push!(h, 1 - ωc / π)
    # Mirror coefficients
    append!(h, reverse(h[1:m]))
end

"    lowpass_coefficients(n; nyquist_denom::Int = 2)"
function lowpass_coefficients(n; nyquist_denom::Int = 2)
    isodd(n) && n > 1 || throw(ArgumentError("n"))
    nyquist_denom > 1 || throw(ArgumentError("nyquist_denom"))
    # Even number of coefficents on one side of center
    m = div(n,  2)
    # Cut-off frequency
    ωc = π / nyquist_denom
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
    freqs = (0:(n - 1)) .* 2π / (n - 1)
    collect(0.42 .- 0.5 * cos.( freqs) .+ 0.08 * cos.(2 * freqs))
end

"""
    fir_lp_coefficients(n; nyquist_denom::Int = 10)

n is window length, an odd number
nyquist_denom, here ν,  is used in calculating cutoff frequency, equivalent to
cutoff wave length when samples are spaced evenly. By cutoff, we mean -6dB attenuation.

    ωc = 2π / ν

In terms of critical wave length:

    λc = 2ν 
"""
function fir_lp_coefficients(n; nyquist_denom::Int = 10)
    c = lowpass_coefficients(n; nyquist_denom) .* blackman_coefficients(n)
    # In order to have zero constant-signal attenuation regardless of length,
    # we scale the coeffients to sum one.
    # Also make this an OffsetArray suitable for kernel use
    centered(c ./ sum(c))
end

"""
    fir_hp_coefficients(n;  nyquist_denom::Int = 10)

n is window length, an odd number
nyquist_denom, here ν,  is used in calculating cutoff frequency, equivalent to
cutoff wave length when samples are spaced evenly. By cutoff, we mean -6dB attenuation.

    ωc = 2π / ν

In terms of critical wave length:

    λc = 2ν 
"""
function fir_hp_coefficients(n;  nyquist_denom::Int = 10)
    c = highpass_coefficients(n; nyquist_denom) .* blackman_coefficients(n)
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
    smooth_surface_fir(z; w = 101, nyquist_denom::Int = 10)

Apply a FIR low-pass filter with Blackman window. The intention is to avoid 
phase-distortion while working with something familiar. A Gaussian filter
likely would do the job "better" (with no overshoot), but the parameters are 
less intuitive for me to quantify. We use the Gaussian 'low-pass' filter for 
general blurring.

w is window length, an odd number.
nyquist_denom
"""
function smooth_surface_fir(z; w = 101, nyquist_denom::Int = 10)
    # Coefficients, including a Blackman window.
    c = Float32.(fir_lp_coefficients(w; nyquist_denom))
    # Elevations smoothed by lp-filter
    imfilter(z, (c, permutedims(c)), FIRTiled())
end
smooth_surface_fir(g::GeoArray; w = 101, nyquist_denom::Int = 10) = smooth_surface_fir(permutedims(g.A[:, :, 1]); w, nyquist_denom)


"""
    roughness_of_surface_fir(z;  w = 49, nyquist_denom::Int = 10)

Apply a FIR high-pass filter with Blackman window. The intention is to avoid 
phase-distortion and to work with something familiar.

w is window length, an odd number.
"""
function roughness_of_surface_fir(z;  w = 49, nyquist_denom::Int = 10)
    # Coefficients, including a Blackman window.
    c = Float32.(fir_hp_coefficients(w; nyquist_denom))
    # Elevations where the smooth variation is removed and
    # the bumps remain.
    imfilter(z, (c, transpose(c)), FIRTiled())
end


"""
    bumpiness(z; w = 49, amp_cut⁻ = 1.5, amp_cut⁺ =  4.6477094f0, z_max = 680f0)
    ---> Matrix{Gray{Float32}} (shape and size like z)

Returns positive values defined by a high-pass filter with length w, Blackman window.
Values defined as non-interesting by the other keywords are set to zero. 

After smoothing, would represent an amplitude of elevation variation inside the forest.
"""
function bumpiness(z; w = 75, amp_cut⁻ = 0.5, amp_cut⁺ =  3f0, z_max = 680f0, nyquist_denom::Int = 15)
    isodd(w) || throw(ArgumentError("w must be odd, not $w"))
    # The local bumps (+ and -)
    r = roughness_of_surface_fir(z; w, nyquist_denom)
    r1 = r .* (z .< z_max)
    # Take absolute value, and drop the large values, which are most likely artifacts
    # or houses. 
    Gray{Float32}.(map(r1) do ρ
        bump = abs(ρ)
        bump < amp_cut⁻ ? 0f0 : bump > amp_cut⁺ ? 0f0 : bump
    end)
end




"""
    bumpy_patch(z, source_indices; cell_count_min_m² = 9062, w = 75, amp_cut⁻ = 0.5, amp_cut⁺ =  3f0, z_max = 680f0, nyquist_denom::Int = 15)
    bumpy_patch(g::GeoArray, source_indices::CartesianIndices)
    ---> typeof(z)

Used by 'ridges' and 'contours'.

# Arguments

z                 Elevations. The other default parameters assumes a grid spacing of 1 utm meter.
cell_count_min    Remove smaller forests
w                 Window length for high-pass filter
amp_cut⁻          Drop filtered values below
amp_cut⁺          Drop filtered values above (often steep ridges or artifacts)
z_max             Drop forest above this elevation
"""
function bumpy_patch(z, source_indices; cell_count_min_m² = 9062, w = 75, amp_cut⁻ = 0.5, amp_cut⁺ =  3f0, z_max = 680f0, nyquist_denom::Int = 15)
    # Isolate short-wavelength variation (mostly forests)
    b = bumpiness(z; w, amp_cut⁻, z_max, amp_cut⁺, nyquist_denom)
    # Forests are interspersed with high and low bumps. Fill the
    # spaces between trees by cutting off the tree tops.
    #
    # After smoothing, we don't need the full resolution anymore:
    c = smooth_surface_fir(b; nyquist_denom = 4, w)[source_indices]
    #
    # Threshold to black-and white. Isolated islands inside bumpy regions are OK.
    # Must be a Float for segmentation...
    d = map(x -> Gray{Float32}(x > 0.01), c)
    # Over-segmented image. The exact parameter value doesn't matter, we currently have only 1 and 0 values.
    g = fast_scanning(d, 0.2)
    # Merge small segments into larger.
    cell2utm = source_indices.indices[1].step
    cell_count_min = div(cell_count_min_m²,  cell2utm^2)
    h = prune_iterations(g, cell_count_min)
    # Change to binary matrix image
    map(labels_map(h)) do j
        amp_cut⁻ < segment_mean(h, j) < amp_cut⁺ ? 1f0 : 0f0
    end
end
# TODO: Use 'elevation_...' function in caller
bumpy_patch(g::GeoArray, source_indices::CartesianIndices) = bumpy_patch(transpose(g.A[:, :, 1]), source_indices)


function prune_iterations(g::SegmentedImage, cell_count_min; iterations = 3)
    h = g
    i = 0
    while i < iterations
        f_to_be_removed = j -> (segment_pixel_count(h, j) < cell_count_min)
        f_diff = (rem_label, neigh_label) -> -segment_mean(h, neigh_label)
        h = prune_segments(h, f_to_be_removed, f_diff)
        i += 1
    end
    h
end