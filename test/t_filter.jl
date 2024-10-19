using Test
using BitmapMaps: fir_hp_coefficients, fir_lp_coefficients, centered, imfilter
using BitmapMaps: highpass_coefficients, lowpass_coefficients, blackman_coefficients
using BitmapMaps: mapwindow, Gray, mark_at!, line!, colorview, display_if_vscode
import BitmapMaps

# Function to calculate the RMS value of a signal
rms(signal) = sqrt(sum(signal .^ 2) / length(signal))
function attenuation_dB(input_signal, output_signal)
    rms_input = rms(input_signal)
    rms_output = rms(output_signal)
    20 * log10(rms_output / rms_input)
end
# Generate signal
"""
    vsin(;n, w, l)

Sinus, n cycles per filter length.
Sum zero for n window lengths in sample length.
Sum +/- 1 for half cycles in sample length.

# Arguments

- n is cycles per window length
- w is window length
- l is sample length
"""
vsin(;n, w, l) = [sin(2π * (n * (i - 1) / w)) for i in 1:l]
@test abs(sum(vsin(n = 1, w = 4, l = 4))) < 0.00001
@test abs(sum(vsin(n = 1, w = 4, l = 6)) - 1) < 0.00001
@test abs(sum(vsin(n = 1, w = 5, l = 50))) < 0.00001

function attenuation_dB(λ::Int; l = 120, f = input -> imfilter(input, fir_hp_coefficients(9)), plt = true)
    n = l / λ
    input_signal = vsin(;n, w = l, l)
    @assert length(input_signal) == l
    output_signal = f(input_signal)
    if plt && @isdefined UnicodePlots
        p = lineplot(input_signal, width = 120, xlim = (1, l))
        lineplot!(p, output_signal)
        println(p)
    end
    attenuation_dB(input_signal, output_signal)
end

function attenuation_dB(λ;  l = 2000, f = input -> imfilter(input, fir_hp_coefficients(9)), plt = true, p = nothing)
    v = attenuation_dB.(λ;  l, f, plt = false)
    if plt && @isdefined UnicodePlots
        kwds = (width = 120, xlabel = "Wave length λ", ylabel = "dB",
            xlim = (1, maximum(λ)), 
            ylim = (floor(minimum(v)), ceil(maximum(v))))
        if isnothing(p)
            p = lineplot(λ, v; kwds...)
        end
        p = lineplot!(p, λ, v)
        hline!(p, [1, 0], [maximum(λ), 0], color=:cyan)
        hline!(p, [1, -6], [maximum(λ), -6], color=:cyan)
    else 
        p = "UnicodePlots not loaded"
    end
    # print a table
    printstyled(lpad("λ Wave length [samples]", 30), lpad("Attenuation [dB]", 30), "\n", color = :green)
    for (len, att) in zip(λ, v)
        if -6 < att < -3
            color = :white
        elseif att < -6
            color = :green 
        else
            color = :light_black
        end
        printstyled(lpad(string(len), 30), lpad(string(round(att, digits = 1)), 30), "\n"; color)
    end
    p
end    
function attenuation_dB_hp_table(; l = 2000, w = 49)
    nyq_denom_range = 2:8
    λ = 1:30
    M = zeros(Float64, length(λ), length(nyq_denom_range))
    for (i, nyquist_denom) in enumerate(nyq_denom_range)
        f = input -> imfilter(input, fir_hp_coefficients(w; nyquist_denom))
        v = attenuation_dB.(λ;  l, f, plt = false)
        M[:, i] = v
    end
    printstyled(lpad("λ Wave length [samples]", 30), lpad("Attenuation [dB] @ nyquist_denom", 50), "\n", color = :white)
    printstyled(lpad(" ", 30), color = :white)
    for n in nyq_denom_range
        printstyled(lpad(string(n), 15), color = :white)
    end
    println()
    for (i, len) in enumerate(λ)
        printstyled(lpad(string(len), 30); color = :normal)
        for (j, nyquist_denom) in enumerate(nyq_denom_range)
            att = M[i, j]
            if -6 < att < -3
                color = :white
            elseif att < -6
                color = :green 
            else
                color = :light_black
            end
            printstyled(lpad(string(round(att, digits = 1)), 15); color)
        end
        println()
    end
end



attenuation_dB_hp_table()
attenuation_dB(1:10)
f = input -> imfilter(input, fir_hp_coefficients(99; nyquist_denom = 4) )
p = attenuation_dB(1:10; f)
f = input -> imfilter(input, fir_hp_coefficients(99; nyquist_denom = 5) )
p = attenuation_dB(1:15; f, p)


f2 = input -> imfilter(input, fir_hp_coefficients(99; nyquist_denom = 2) )
f3 = input -> imfilter(input, fir_hp_coefficients(99; nyquist_denom = 3) )
f4 = input -> imfilter(input, fir_hp_coefficients(99; nyquist_denom = 4) )
f5 = input -> imfilter(input, fir_hp_coefficients(99; nyquist_denom = 5) )
f6 = input -> imfilter(input, fir_hp_coefficients(99; nyquist_denom = 6) )
f7 = input -> imfilter(input, fir_hp_coefficients(99; nyquist_denom = 7) )
f8 = input -> imfilter(input, fir_hp_coefficients(99; nyquist_denom = 8) )

p = attenuation_dB(1:25; f = f2)
attenuation_dB(1:25; f = f3, p)
attenuation_dB(1:25; f = f4, p)
attenuation_dB(1:25; f = f5, p)
attenuation_dB(1:25; f = f6, p)
attenuation_dB(1:25; f = f7, p)
attenuation_dB(1:25; f = f8, p)



function attenuation_dB_lp_table(; l = 2000, w = 49)
    nyq_denom_range = 2:8
    λ = 1:30
    M = zeros(Float64, length(λ), length(nyq_denom_range))
    for (i, nyquist_denom) in enumerate(nyq_denom_range)
        f = input -> imfilter(input, fir_lp_coefficients(w; nyquist_denom))
        v = attenuation_dB.(λ;  l, f, plt = false)
        M[:, i] = v
    end
    printstyled(lpad("λ Wave length [samples]", 30), lpad("Attenuation [dB] @ nyquist_denom", 50), "\n", color = :white)
    printstyled(lpad(" ", 30), color = :white)
    for n in nyq_denom_range
        printstyled(lpad(string(n), 15), color = :white)
    end
    println()
    for (i, len) in enumerate(λ)
        printstyled(lpad(string(len), 30); color = :normal)
        for (j, nyquist_denom) in enumerate(nyq_denom_range)
            att = M[i, j]
            if -6 < att < -3
                color = :white
            elseif att < -6
                color = :green 
            else
                color = :light_black
            end
            printstyled(lpad(string(round(att, digits = 1)), 15); color)
        end
        println()
    end
end

attenuation_dB_lp_table()
attenuation_dB_lp_table(; w = 15)