# This file is a part of RadiationSpectra.jl, licensed under the MIT License (MIT).
using Test

using Distributions
using StatsBase
using IntervalSets
using ValueShapes

# using BAT
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

    @testset "Custom Spectrum Density" begin
        struct CustomSpectrumDensity{T}  <: RadiationSpectra.UvSpectrumDensity{T}
            A::T
            σ::T
            μ::T
            offset::T
        end
        RadiationSpectra.evaluate(d::CustomSpectrumDensity, x) = d.A * pdf(Normal(d.μ, d.σ), x) + d.offset
        # Constructor
        function CustomSpectrumDensity(nt::NamedTuple{(:A, :μ, :σ, :offset)}) 
            T = promote_type(typeof.(values(nt))...)
            CustomSpectrumDensity(T(nt.A), T(nt.σ), T(nt.μ), T(nt.offset))
        end
        p0 = (A = 0.9*N, μ = 0.0, σ = 2.0, offset = 1.0)
        lower_bounds = (A = 0.0, μ = -5, σ = 0.01, offset = 0.0)
        upper_bounds = (A = 10.0*N, μ = 5, σ = 200.0, offset = 500.0)

        fitted_dist, backend_result = fit(CustomSpectrumDensity, h, p0, lower_bounds, upper_bounds)
        @test isapprox(true_pars[1], fitted_dist.A, rtol = 1e-2)
        @test isapprox(true_pars[2], fitted_dist.μ, rtol = 1e-2)
        @test isapprox(true_pars[3], fitted_dist.σ, rtol = 1e-2)
        @test isapprox(0, fitted_dist.offset, atol = N / 1000)
        
        # Fix one parameter in the fit:
        p0 = (A = 0.9*N, μ = 0.0, σ = 2.0, offset = 0)
        shape = NamedTupleShape(
            A = ScalarShape{Real}(),
            μ = ScalarShape{Real}(),
            σ = ScalarShape{Real}(),
            offset = ConstValueShape{Real}(p0.offset),
        )
        fitted_dist, backend_result = fit(CustomSpectrumDensity, h, p0, lower_bounds, upper_bounds, shape)
        @test isapprox(true_pars[1], fitted_dist.A, rtol = 1e-2)
        @test isapprox(true_pars[2], fitted_dist.μ, rtol = 1e-2)
        @test isapprox(true_pars[3], fitted_dist.σ, rtol = 1e-2)
        @test fitted_dist.offset == p0.offset
    end

    include("autocalibration.jl")

    # @testset "BAT Backend" begin
    #     # ENV["JULIA_DEBUG"] = "BAT"
    #     fitted_dist, bat_samples_result = RadiationSpectra.fit(RadiationSpectra.NormalPeakUvD, float(h), backend = :BAT, strict = false)
    #     @test isapprox(true_pars[1], fitted_dist.A, rtol = 1e-2)
    #     @test isapprox(true_pars[2], fitted_dist.UvNormal.μ, rtol = 1e-2)
    #     @test isapprox(true_pars[3], fitted_dist.UvNormal.σ, rtol = 1e-2)
    # end
end # testset