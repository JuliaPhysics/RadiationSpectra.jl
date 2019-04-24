function _leastsquare(ydata::Vector{T}, ymodel::Vector{T}, weights::Vector{T})::T where {T <: AbstractFloat}
    ls::T = 0
    @inbounds for i in eachindex(ydata)
        ls += (ydata[i] - ymodel[i])^2 * weights[i]
    end
    return ls
end

function lsqfit!(fit::FitFunction{T, 1, NP}, xdata::Vector{T}, ydata::Vector{T}; weights::Vector{T} = ones(T, length(xdata))) where {T <: AbstractFloat, NP}
    function lsq(params::Vector{T})::T
        ymodel::Vector{T} = fit.model(xdata, params)::Vector{T}
        ls::T = _leastsquare( ydata, ymodel, weights )
        return ls
    end 
    optim_result = try
        optimize( lsq, fit.initial_parameters) 
    catch err
        optimize( lsq, fit.initial_parameters, BFGS())
    end
    _set_fitted_parameters!(fit, optim_result.minimizer)
end

function lsqfit!(fit::FitFunction{T, 1, NP}, xdata::Vector{T}, ydata::Vector{T}, xerr::Vector{T}, yerr::Vector{T}) where {T <: AbstractFloat, NP}
    weights::Vector{T} = @. 1 / sqrt(xerr^2 + yerr^2)
    return lsqfit!(fit, xdata, ydata, weights = weights)
end

function lsqfit!(fit::FitFunction{T, 1, NP}, xdata::Vector{T}, ydata::Vector{T}, err::Vector{T}) where {T <: AbstractFloat, NP}
    weights::Vector{T} = inv.(err)
    return lsqfit!(fit, xdata, ydata, weights = weights)
end


"""
    lsqfit!(fit::FitFunction{T, 1, NP}, h::Histogram)::Nothing where {T <: AbstractFloat, NP}

Performs a least square fit with the model `fit.model` and the initial parameters `fit.initial_parameters`
on the histogram `h` in the range `fit.fitranges[1]`. The determined parameters are stored in `fit.fitted_parameters`.
"""
function lsqfit!(fit::FitFunction{T, 1, NP}, h::Histogram)::Nothing where {T <: AbstractFloat, NP}
    first_bin::Int = StatsBase.binindex(h, first(fit.fitranges[1]))
    last_bin::Int  = StatsBase.binindex(h, last(fit.fitranges[1]))
    if first_bin < 1 first_bin = 1 end
    if (last_bin > length(h.weights)) last_bin = length(h.weights) end
    bin_centers::Vector{T} = StatsBase.midpoints(h.edges[1])[first_bin:last_bin]
    counts::Vector{T} = h.weights[first_bin:last_bin]
    weights::Vector{T} = sqrt.(counts) # Poisson distributed
    weights = [ w != 0 ? w : 1.  for w in weights] 

    lsqfit!(fit, bin_centers, counts, weights)
end


