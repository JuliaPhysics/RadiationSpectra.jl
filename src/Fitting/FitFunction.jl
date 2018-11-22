"""
    AbstractFitFunction{T, N}

Abstract type for an `N`-dimensional fit of eltype `T`.
"""
abstract type AbstractFitFunction{T, N} end

"""
    FitFunction{T} <: AbstractFitFunction{T, 1}

Fields:
- `model::Function`: Function of the fit model.
- `fitrange::Tuple{Union{Missing, T}, Union{Missing, T}}`: Range on which the fit is performed.
- `parameters::Union{Vector{Missing}, Vector{T}}`: Fitted parameters.
- `uncertainties::Union{Vector{Missing}, Vector{T}}`: Uncertainties of the fitted parameters.
- `initial_parameters::Union{Vector{Missing}, Vector{T}}`: Initial parameters.
- `confidence_level::Union{Missing, T}`: Confidence level, used to determine the uncertainties.

# Plotting:
A plot recipes exists for this struct: `plot(fit::FitFunction{T}; npoints=501)`

Plots the model function over the fitrange with 501-points. 
"""
mutable struct FitFunction{T} <: AbstractFitFunction{T, 1}
    model::Function
    fitrange::Tuple{Union{Missing, T}, Union{Missing, T}}
    parameters::Union{Vector{Missing}, Vector{T}}
    uncertainties::Union{Vector{Missing}, Vector{T}}
    initial_parameters::Union{Vector{Missing}, Vector{T}}
    confidence_level::Union{Missing, T}

    function FitFunction(model::Function, fitrange::AbstractRange, parameters::Vector{T}, uncertainties::Vector{T}, initial_parameters::Vector{T}, confidence_level::T) where {T <: Real}
        return new{T}(model, fitrange, parameters, uncertainties, initial_parameters, confidence_level)
    end
    function FitFunction( model::Function )
        return new{Float64}( model, (0,1), zeros(Float64, 0), zeros(Float64, 0), zeros(Float64, 0), 0.68)
    end
end


function print(io::IO, f::FitFunction)
    print(io, "Model: ", f.model, " - fit range: ", f.fitrange, " - parameters: ", f.parameters, " - uncertainties: ", f.uncertainties, " - initial_parameters: ", f.initial_parameters, " - confidence_level: ", f.confidence_level )
end
function println(io::IO, f::FitFunction)
    println(io, "Model: ", f.model)
    println(io, "fit range: ", f.fitrange)
    println(io, "parameters: ", f.parameters)
    println(io, "initial_parameters: ", f.initial_parameters)
    println(io, "uncertainties: ", f.uncertainties)
    println(io, "confidence_level: ", f.confidence_level)
end

function show(io::IO, f::FitFunction) 
    println(io, f)
end
function show(io::IO, ::MIME"text/plain", f::FitFunction) 
    show(io, f)
end

