# This file is a part of RadiationSpectra.jl, licensed under the MIT License (MIT).
using Test

using Distributions
using StatsBase
using IntervalSets
using ValueShapes

using BAT
using RadiationSpectra

@testset "Package RadiationSpectra" begin
    T = Float64
    N = 10^6
    true_pars = T[N, 0.5, 1]
    h = fit(Histogram, rand(Normal(true_pars[2:3]...), N), nbins = 400)
     
    @testset "Fit - Optimizer" begin 
        fitted_dist, opt_result_nt = RadiationSpectra.fit(RadiationSpectra.NormalPeakUvD, h) 
        @test isapprox(true_pars[1], fitted_dist.A, rtol = 1e-2)
        @test isapprox(true_pars[2], fitted_dist.UvNormal.μ, rtol = 1e-2)
        @test isapprox(true_pars[3], fitted_dist.UvNormal.σ, rtol = 1e-2)
        @test isapprox(0, fitted_dist.bgll, rtol = true_pars[1] * 1e-2)
        @test isapprox(0, fitted_dist.bglr, rtol = true_pars[1] * 1e-2)
    end

    @testset "Auto Calibration" begin
        h_uncal = RadiationSpectra.get_example_spectrum()
        
        gamma_lines = T[609.312, 911.204, 1120.287, 1460.830, 1764.494, 2614.533]
        h_cal, h_deconv, peakPositions, threshold, c, c_precal = RadiationSpectra.calibrate_spectrum(h_uncal, gamma_lines, min_n_peaks = 30, σ = 2.0)
        
        # @info "Calibration constant: c = $(c)"
        c_true = 0.011259105696794384
        @testset "Auto calibration" begin
            @test abs(c - c_true) / c_true < 0.001
        end
    end
    
    @testset "BAT Backend" begin
        # ENV["JULIA_DEBUG"] = "BAT"
        fitted_dist, bat_samples_result = RadiationSpectra.fit(RadiationSpectra.NormalPeakUvD, float(h), backend = :BAT, strict = false);
        @test isapprox(true_pars[1], fitted_dist.A, rtol = 1e-2)
        @test isapprox(true_pars[2], fitted_dist.UvNormal.μ, rtol = 1e-2)
        @test isapprox(true_pars[3], fitted_dist.UvNormal.σ, rtol = 1e-2) 
    end
end # testset

