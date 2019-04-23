# This file is a part of RadiationSpectra.jl, licensed under the MIT License (MIT).
using RadiationSpectra

using Test
@testset "Package RadiationSpectra" begin
    T = Float64

    h_uncal = RadiationSpectra.get_example_spectrum()

    gamma_lines = T[609.312, 668, 785., 911.204, 1120.287, 1460.830, 1764.494, 2614.533]
    h_cal, h_deconv, peakPositions, threshold, c, c_precal = RadiationSpectra.calibrate_spectrum(h_uncal, gamma_lines, min_n_peaks = 30, σ = 2.0 )

    @info "Calibration constant: c = $(c)"
    c_true = 0.011269787869370062
    @testset "Auto calibration" begin
        @test abs(c - c_true) / c_true < 0.025
    end

    ff = FitFunction(T, :GaussPlusLinearBackground)

    set_fitranges!(ff, ((1461 - 20, 1461 + 20),))
    set_initial_parameters!(ff, T[100, 1, 1461, 0, 0])
    lsqfit!(ff, h_cal)

    fitted_pars = get_fitted_parameters(ff)
    @show fitted_pars
    @testset "LSQ Fit" begin
        @test abs(fitted_pars[:μ] - 1460.830) <= 3
    end
    
    set_initial_parameters!(ff, fitted_pars)
    llhfit!(ff, h_cal)
    
    fitted_pars = get_fitted_parameters(ff)
    @show fitted_pars
    @testset "LLH Fit" begin
        @test abs(fitted_pars[:μ] - 1460.830) <= 3
    end


end # testset
