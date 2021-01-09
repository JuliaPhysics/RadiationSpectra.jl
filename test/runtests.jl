# This file is a part of RadiationSpectra.jl, licensed under the MIT License (MIT).
using Test

using Distributions
using StatsBase
 
using RadiationSpectra

@testset "Package RadiationSpectra" begin
    T = Float64
    N = 10^6
    true_pars = T[N, 0.5, 1]
    h = fit(Histogram, rand(Normal(true_pars[2:3]...), N), nbins = 400)
     
    @testset "rsfit - Parameter as NamedTuple" begin 
        fitted_pars_nt, opt_result_nt = RadiationSpectra.rsfit(RadiationSpectra.NormalPeakUvD, h) 
        @test isapprox(true_pars[1], fitted_pars_nt.A, rtol = 1e-2)
        @test isapprox(true_pars[2], fitted_pars_nt.μ, rtol = 1e-2)
        @test isapprox(true_pars[3], fitted_pars_nt.σ, rtol = 1e-2)
        @test isapprox(0, fitted_pars_nt.bgll, rtol = true_pars[1] * 1e-2)
        @test isapprox(0, fitted_pars_nt.bglr, rtol = true_pars[1] * 1e-2)
    end

    @testset "Auto Calibration" begin
        h_uncal = RadiationSpectra.get_example_spectrum()
        
        gamma_lines = T[609.312, 911.204, 1120.287, 1460.830, 1764.494, 2614.533]
        h_cal, h_deconv, peakPositions, threshold, c, c_precal = RadiationSpectra.calibrate_spectrum(h_uncal, gamma_lines, min_n_peaks = 30, σ = 2.0 )
        
        @info "Calibration constant: c = $(c)"
        c_true = 0.011259105696794384
        @testset "Auto calibration" begin
            @test isapprox(c, c_true, rtol = 1e-3) #abs(c - c_true) / c_true < 0.001
        end
    end

    
    @testset "BAT Backend" begin
        # ENV["JULIA_DEBUG"] = "BAT"
        # using BAT
        # bat_samples = RadiationSpectra.rsbatfit(RadiationSpectra.NormalPeakUvD, h) 

        # using ForwardDiff
        # using ValueShapes

        # T = Float64
        # bounds = RadiationSpectra.initial_parameter_guess(h, RadiationSpectra.NormalPeakUvD{T})[end]
        # hp = RadiationSpectra.BATHistLLHPrecalulations(RadiationSpectra.HistLLHPrecalulations(h), RadiationSpectra.NormalPeakUvD)
        # prior = BAT.DistributionDensity(BAT.NamedTupleDist( bounds ) )
        # globalparshape = varshape(prior)
        # x = bat_initval(prior).result
    
        # BAT.eval_logval_unchecked(hp, x[])
    
        # @time BAT.eval_logval_unchecked(hp, x[])
    
        # X = unshaped(x, globalparshape);
        # f = let hp=hp, shape=globalparshape
        #     X -> BAT.eval_logval_unchecked(hp, shape(X)[])
        # end
    
        # ForwardDiff.gradient(f, X)
    
        # BAT.eval_gradlogval(hp, x)
    
        
    end


end # testset

