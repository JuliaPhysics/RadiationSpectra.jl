# This file is a part of RadiationSpectra.jl, licensed under the MIT License (MIT).
using RadiationSpectra

import Test
Test.@testset "Package RadiationSpectra" begin
    T = Float64

    h_uncal = RadiationSpectra.get_example_spectrum()

    gamma_lines = T[609.312, 668, 785., 911.204, 1120.287, 1460.830, 1764.494, 2614.533]
    h_cal, h_deconv, peakPositions, threshold, c, c_precal = RadiationSpectra.calibrate_spectrum(h_uncal, gamma_lines, min_n_peaks = 30, σ = 2.0 )

    @info "Calibration constant: c = $(c)"
    c_true = 0.011269787869370062
    Test.@test abs(c - c_true) / c_true < 0.025

    ff = FitFunction(T, :GaussPlusLinearBackground)

    set_fitranges!(ff, ((1461 - 20, 1461 + 20),))
    set_initial_parameters!(ff, T[100, 1, 1461, 0, 0])
    lsqfit!(ff, h_cal)
    Test.@test abs(get_fitted_parameters(ff)[:μ] - 1460.830) <= 3

end # testset
