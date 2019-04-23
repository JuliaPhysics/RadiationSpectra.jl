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

    ff_lsq = FitFunction(T, :GaussPlusLinearBackground)

    set_fitranges!(ff_lsq, ((1461 - 20, 1461 + 20),))
    set_initial_parameters!(ff_lsq, T[10000, 1.4, 1461, 20, 0])
    lsqfit!(ff_lsq, h_cal)
    
    fitted_pars = collect(get_fitted_parameters(ff_lsq))
    @testset "LSQ Fit" begin
        @test abs(fitted_pars[3] - 1460.830) <= 3
    end

    ff_llh = FitFunction(T, :GaussPlusLinearBackground)
    set_fitranges!(ff_llh, ((1461 - 20, 1461 + 20),))
    fitted_pars[1] = 0.9 * fitted_pars[1]

    set_initial_parameters!(ff_llh, fitted_pars)
    llhfit!(ff_llh, h_cal)
    
    fitted_pars = get_fitted_parameters(ff_llh)
    @show ff_llh
    @testset "LLH Fit" begin
        @test abs(fitted_pars[3] - 1460.830) <= 3
    end


    @testset "General model functions" begin
        ff = FitFunction(T, :Gauss)
        RadiationSpectra._set_fitted_parameters!(ff, [1, 1, 0])
        ff.model(T[0, 1], collect(get_fitted_parameters(ff)))
        println(ff)
        @test round.(ff.model(T[0, 1], collect(get_fitted_parameters(ff))), digits = 5) == round.(T[0.3989422804014327, 0.24197072451914337], digits = 5)
        
        ff = FitFunction(:PDF_Gauss)
        ff = FitFunction(T, :PDF_Gauss)
        RadiationSpectra._set_fitted_parameters!(ff, [1, 0])
        ff.model(T[0, 1], collect(get_fitted_parameters(ff)))
        println(ff)
        @test round.(ff.model(T[0, 1], collect(get_fitted_parameters(ff))), digits = 5) == round.(T[0.3989422804014327, 0.24197072451914337], digits = 5)
        
        ff = FitFunction(T, :GaussPlusConstantBackground)
        RadiationSpectra._set_fitted_parameters!(ff, [1, 1, 0, 0])
        ff.model(T[0, 1], collect(get_fitted_parameters(ff)))
        println(ff)
        @test round.(ff.model(T[0, 1], collect(get_fitted_parameters(ff))), digits = 5) == round.(T[0.3989422804014327, 0.24197072451914337], digits = 5)
        
        ff = FitFunction(T, :GaussPlusLinearBackground)
        RadiationSpectra._set_fitted_parameters!(ff, [1, 1, 0, 0, 0])
        ff.model(T[0, 1], collect(get_fitted_parameters(ff)))
        println(ff)
        @test round.(ff.model(T[0, 1], collect(get_fitted_parameters(ff))), digits = 5) == round.(T[0.3989422804014327, 0.24197072451914337], digits = 5)
    end
    

end # testset
