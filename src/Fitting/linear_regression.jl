@fastmath function linear_regression(x::Vector{<:Real}, y::Vector{<:Real})::Vector{T} # Substitutes linear fit --> much faster
    @assert length(x) == length(y) "x and y must have the same length."
    T=Float64
    x_mean::T = mean(x)
    y_mean::T = mean(y)
    num::T = 0.0
    nom::T = 0.0
    for i in eachindex(x)
        x_res = (x[i] - x_mean)
        num += x_res * (y[i] - y_mean)
        nom += x_res*x_res
    end
    slope::T = num/nom
    offset::T = y_mean - slope * x_mean
    return [offset, slope]
end
@fastmath function linear_regression(x::Vector{T}, y::Vector{T})::Vector{T} where {T <: AbstractFloat} # Substitutes linear fit --> much faster
    @assert length(x) == length(y) "x and y must have the same length."
    x_mean::T = mean(x)
    y_mean::T = mean(y)
    num::T = 0.0
    nom::T = 0.0
    @inbounds for i in eachindex(x)
        x_res = (x[i] - x_mean)
        num += x_res * (y[i] - y_mean)
        nom += x_res * x_res
    end
    slope::T = num / nom
    offset::T = y_mean - slope * x_mean
    return [offset, slope]
end


function linear_regression!(fit::FitFunction{T, 1, NP}, xdata::Vector{T}, ydata::Vector{T}; weights::Vector{T} = ones(T, length(xdata)), kwargs...) where {T <: AbstractFloat, NP}
    @assert fit.model == RadiationSpectra.Linear
    fit.fitranges == ([-Inf, Inf],) ? set_fitranges!(fit, ((xdata[1], xdata[end]), )) : nothing
    param = linear_regression(xdata, ydata)
    _set_fitted_parameters!(fit, param)
    _set_residuals!(fit, xdata, ydata)
    _set_Χ²!(fit, xdata, fit.residuals)
    fit
end
