using Test
using BitmapMaps: fir_hp_coefficients, fir_lp_coefficients, centered, imfilter
using BitmapMaps: highpass_coefficients, lowpass_coefficients, blackman_coefficients
using BitmapMaps: mapwindow
import BitmapMaps
@testset "Filter Coefficient Tests" begin
    @test length(fir_hp_coefficients(9)) == 9
    @test length(fir_lp_coefficients(9)) == 9
    @test fir_hp_coefficients(3) ≈ centered([8.617992180508635e-19, 0.9374999999999999, 8.617992180508635e-19])
    @test fir_lp_coefficients(3) ≈ centered([-1.3788787488813818e-17, 1.0, -1.3788787488813818e-17])
    @test_throws AssertionError fir_hp_coefficients(2)
    @test_throws AssertionError fir_lp_coefficients(1)
    @test abs(sum(fir_hp_coefficients(99))) < 0.001
    @test abs(sum(fir_hp_coefficients(9))) < 0.8
    @test sum(fir_lp_coefficients(99)) ≈ 1
    @test sum(fir_lp_coefficients(9)) ≈ 1
end


# Constant signal
@test sum(imfilter(ones(1000), fir_hp_coefficients(9))) < 800
@test sum(imfilter(ones(1000), fir_hp_coefficients(27))) < 390
@test sum(imfilter(ones(1000), fir_hp_coefficients(81))) < 2
@test sum(imfilter(ones(1000), fir_lp_coefficients(9))) ≈ 1000
@test sum(imfilter(ones(1000), fir_lp_coefficients(27))) ≈ 1000
@test sum(imfilter(ones(1000), fir_lp_coefficients(81))) ≈ 1000
# Ramp signal
rampe = range(0, 1, length = 2000)
@test sum(rampe) == 1000
@test sum(imfilter(rampe, fir_hp_coefficients(9))) < 800
@test sum(imfilter(rampe, fir_hp_coefficients(27))) < 390
@test sum(imfilter(rampe, fir_hp_coefficients(81))) < 2
@test sum(imfilter(rampe, fir_lp_coefficients(9))) == 1000
@test sum(imfilter(rampe, fir_lp_coefficients(27))) ≈ 1000 # There's a small frequncy contribution at flat end
@test sum(imfilter(rampe, fir_lp_coefficients(81))) ≈ 1000
# Function to calculate the RMS value of a signal
rms(signal) = sqrt(sum(signal .^ 2) / length(signal))
function attenuation_dB(input_signal, output_signal)
    rms_input = rms(input_signal)
    rms_output = rms(output_signal)
    20 * log10(rms_output / rms_input)
end
# Sinus, n cycles per filter length
# Sum zero for n window lengths in sample length
# Sum +/- 1 for half cycles in sample length
# n is cycles per window length
# w is window length
# l is sample length
vsin(;n, w, l) = [sin(2π * (n * (i - 1) / w)) for i in 1:l]
@test abs(sum(vsin(n = 1, w = 4, l = 4))) < 0.00001
@test abs(sum(vsin(n = 1, w = 4, l = 6)) - 1) < 0.00001
@test abs(sum(vsin(n = 1, w = 5, l = 50))) < 0.00001
# Half sinus integrated over window length is 2.
n = 0.5; w = 9; l = 9
@test abs(sum(vsin(; n, w, l)) / w - 2 / π) < 0.01
# This half-wave is not much different
@test abs(sum(imfilter(vsin(;n, w, l), fir_lp_coefficients(w))) / w - 2 / π) < 0.02

# Filter attenuation, half of Nyquist frequency
n = 9 / 4; w = 9; l = 900
input = vsin(;n, w, l)
output = imfilter(input, fir_lp_coefficients(w))
@test abs(rms(input) - √(2) / 2) < 0.0000001
@test rms(output) / rms(input) < 0.11
@test attenuation_dB(input, output) < -19.5

# Filter attenuation, full wave over window, 
n = 1; w = 9; l = 900
input = vsin(;n, w, l)
output = imfilter(input, fir_lp_coefficients(w))
@test rms(output) / rms(input) < 0.68
@test attenuation_dB(input, output) < -3.4

# Since our coeffients are symmetric, the 'convolution' vs 'correlation' distiction does not matter:
@test imfilter(input, BitmapMaps.ImageFiltering.reflect(fir_lp_coefficients(w))) == output

# Implement same with 'mapwindow'
function func_lp(; w = 9)
    c = collect(fir_lp_coefficients(w)) # maybe svec.
    mid = w ÷ 2 + 1
    (z) -> begin
        zf = BitmapMaps.conv(z, c)
        zf[mid]
    end
end
flp = func_lp(; w)
output1 = mapwindow(flp, input, (9, ))
@test output == output1
# using BenchmarkTools
#   4.457 μs (10 allocations: 22.08 KiB)
# @btime imfilter(input, $(fir_lp_coefficients(w)))  # Algorith FIRTiled probably even faster
# 116.300 μs (1868 allocations: 292.95 KiB)
# @btime mapwindow(flp, input, (9, ))

# NOTE: Since 'imfilter' is ridiculously advanced and fast, stick to using that even though
# we need to allocate big temporary arrays for the derivatives during calculation of contours.
