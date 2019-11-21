function lsqfit!(fit::FitFunction{T, 1, NP}, xdata::Vector{T}, ydata::Vector{T}; weights::Vector{T} = ones(T, length(xdata)), kwargs...) where {T <: AbstractFloat, NP}
    f = curve_fit(fit.model, xdata, ydata, weights, fit.initial_parameters; kwargs...)
    set_fit_backend_result!(fit, f)
    _set_fitted_parameters!(fit, f.param)
    fit
end

function lsqfit!(fit::FitFunction{T, 1, NP}, xdata::Vector{T}, ydata::Vector{T}, xerr::Vector{T}, yerr::Vector{T}; kwargs...) where {T <: AbstractFloat, NP}
    weights::Vector{T} = @. 1 / (xerr^2 + yerr^2)
    return lsqfit!(fit, xdata, ydata, weights = weights; kwargs...)
end

function lsqfit!(fit::FitFunction{T, 1, NP}, xdata::Vector{T}, ydata::Vector{T}, err::Vector{T}; kwargs...) where {T <: AbstractFloat, NP}
    weights::Vector{T} = inv.(err.^2)
    return lsqfit!(fit, xdata, ydata, weights = weights; kwargs...)
end


"""
    lsqfit!(fit::FitFunction{T, 1, NP}, h::Histogram)::Nothing where {T <: AbstractFloat, NP}

Performs a least square fit with the model `fit.model` and the initial parameters `fit.initial_parameters`
on the histogram `h` in the range `fit.fitranges[1]`. The determined parameters are stored in `fit.fitted_parameters`.
"""
function lsqfit!(fit::FitFunction{T, 1, NP}, h::Histogram) where {T <: AbstractFloat, NP}
    first_bin::Int = StatsBase.binindex(h, first(fit.fitranges[1]))
    last_bin::Int  = StatsBase.binindex(h, last(fit.fitranges[1]))
    if first_bin < 1 first_bin = 1 end
    if (last_bin > length(h.weights)) last_bin = length(h.weights) end
    bin_centers::Vector{T} = StatsBase.midpoints(h.edges[1])[first_bin:last_bin]
    counts::Vector{T} = h.weights[first_bin:last_bin]
    err::Vector{T} = sqrt.(counts) # Poisson distributed
    err = [ w != 0 ? w : 1.  for w in err] 

    lowerbounds = map(b->b.left,  fit.parameter_bounds)
    upperbounds = map(b->b.right, fit.parameter_bounds)

    lsqfit!(fit, bin_centers, counts, err; lower = lowerbounds, upper = upperbounds)
end

_get_standard_deviations(f::FitFunction, fr::LsqFit.LsqFitResult) = LsqFit.stderror(fr)

get_fit_backend_result(fr::LsqFit.LsqFitResult) = fr