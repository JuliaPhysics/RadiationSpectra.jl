# This file is a part of RadiationSpectra.jl, licensed under the MIT License (MIT).

__precompile__(true)

module RadiationSpectra

using DelimitedFiles
using Statistics

using Distributions
using LsqFit
using Optim
using RecipesBase
using SpecialFunctions
using StatsBase

import Base: print, println, show

export peakfinder

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


include("PeakFinder/PeakFinder.jl")

include("Fitting/FitFunction.jl")
include("Fitting/plotrecipes.jl")
include("Fitting/lsqfit.jl")
include("Fitting/llhfit.jl")
include("Fitting/general_model_functions.jl")

include("AutoCalibration/AutoCalibration.jl")


end # module
