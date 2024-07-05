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
    # The sum of coefficients for a long window is zero.
    centered(highpass_coefficients(n) .* blackman_coefficients(n))
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