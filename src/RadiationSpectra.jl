# This file is a part of RadiationSpectra.jl, licensed under the MIT License (MIT).

__precompile__(true)

module RadiationSpectra

using DelimitedFiles
using Statistics

using Distributions
using Optim
using RecipesBase
using StatsBase

import Base: print, println, show

export FitFunction
export get_ndims, get_nparams
export set_fitranges!
export set_parameter_names!
export set_initial_parameters!
export get_fitted_parameters
export get_initial_parameters

export llhfit!, lsqfit!

export peakfinder
export calibrate_spectrum


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
include("Fitting/general_model_functions.jl")

include("Fitting/llhfit.jl")

include("Fitting/lsqfit.jl")

include("AutoCalibration/AutoCalibration.jl")


end # module
