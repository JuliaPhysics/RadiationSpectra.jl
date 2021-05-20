# This file is a part of RadiationSpectra.jl, licensed under the MIT License (MIT).

__precompile__(true)

module RadiationSpectra

    using DelimitedFiles
    using Distributions
    using LinearAlgebra
    using Statistics

    using EmpiricalDistributions
    using IntervalSets
    using NamedTupleTools
    using Optim
    using Random123
    using RecipesBase
    using Requires
    using SpecialFunctions
    using StaticArrays
    using StaticUnivariatePolynomials
    using StatsBase
    using ValueShapes

    import StatsBase: fit

    """
        get_example_spectrum()::Histogram

    Returns an uncalibrated radiation spectrum for testing and demonstrating purpose.
    """
    function get_example_spectrum()::Histogram
        weights::Array{Int, 1} = readdlm(joinpath(@__DIR__, "..", "examples", "data", "hpge-spectrum.txt"))[:]
        h = Histogram(0:100:100 * length(weights), :left)
        h.weights = weights
        return h
    end

    include("general_tools.jl")

    include("ModelDensities.jl")
    
    include("fit.jl")
    include("UnivariateFit.jl")
    
    include("UvModelDensities/NormalPeakUvD.jl")

    include("PeakFinder/PeakFinder.jl")
    include("AutoCalibration/AutoCalibration.jl")

    function __init__()
        @require BAT="c0cd4b16-88b7-57fa-983b-ab80aecada7e" include("BAT.jl")
    end

    include("plot_recipes.jl")

end # module
