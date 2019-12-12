function llhfit(fit::FitFunction{T, 1, NP}, h::Histogram)::Optim.MultivariateOptimizationResults where {T <: AbstractFloat, NP}
    first_bin::Int = StatsBase.binindex(h, first(fit.fitranges[1]))
    last_bin::Int  = StatsBase.binindex(h, last(fit.fitranges[1]))
    if first_bin < 1 first_bin = 1 end
    if (last_bin > length(h.weights)) last_bin = length(h.weights) end
    bin_centers::Vector{T} = StatsBase.midpoints(h.edges[1])[first_bin:last_bin]
    bin_widths::Vector{T} = [StatsBase.binvolume(h, StatsBase.binindex(h, mp)) for mp in bin_centers]
    counts::Vector{T} = h.weights[first_bin:last_bin]

    function log_likelihood(params)
        s = 0
        if in(0, in.(params, fit.parameter_bounds))
            return Inf
        end
        @inbounds for i in eachindex(counts)
            expected_counts = fit.model(bin_centers[i], params)
            if isnan(expected_counts) || expected_counts < 0 
                expected_counts = T(Inf)
            end
            s += -logpdf(Poisson(expected_counts), counts[i])
        end
        return s
    end
    optim_result = Optim.optimize( log_likelihood, fit.initial_parameters )

    uncertainties = try
        sqrt.(diag(inv(ForwardDiff.hessian(log_likelihood, Optim.minimizer(optim_result)))))
    catch e
        zeros(T, length(fit.initial_parameters))
    end
    
    set_fit_backend_result!(fit, (optim_result, uncertainties))
    optim_result
end


"""
    llhfit!(fit::FitFunction{T, 1, NP}, h::Histogram)::Nothing where {T <: AbstractFloat, NP}

Performs a log-likelihood fit of the model function `fit.model` and the initial parameters `fit.initial_parameters`
on the histogram `h` in the range `fit.fitranges[1]`. The determined parameters are stored in `fit.fitted_parameters`.

The likelihood for each individual bin is the Poission distribution.
"""
function llhfit!(fit::FitFunction{T, 1, NP}, h::Histogram) where {T <: AbstractFloat, NP}
    optim_result = llhfit(fit, h)
    _set_fitted_parameters!(fit, optim_result.minimizer)
    fit
end

get_fit_backend_result(fr::Tuple{<:Optim.MultivariateOptimizationResults,<:Any}) = fr[1]

_get_standard_deviations(f::FitFunction, fr::Tuple{<:Optim.MultivariateOptimizationResults,<:Any}) = fr[2] 
