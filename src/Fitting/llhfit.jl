function llhfit(fit::FitFunction, h::Histogram)::Optim.MultivariateOptimizationResults
    first_bin::Int = StatsBase.binindex(h, first(fit.fitrange))
    last_bin::Int  = StatsBase.binindex(h, last(fit.fitrange))
    if first_bin < 1 first_bin = 1 end
    if (last_bin > length(h.weights)) last_bin = length(h.weights) end
    bin_centers::Vector{Float64} = StatsBase.midpoints(h.edges[1])[first_bin:last_bin]
    bin_widths::Vector{Float64} = [StatsBase.binvolume(h, StatsBase.binindex(h, mp)) for mp in bin_centers]
    counts::Vector{Float64} = h.weights[first_bin:last_bin]

    function log_likelihood(params::AbstractVector{Float64})::Float64
        function λ(x, bw) 
            l = fit.model(x, params) #* bw
            return l
        end
        function bin_ll(x, bw, k) 
            try 
                return -logpdf(Poisson(λ(x, bw)), k) 
            catch err 
                return Inf
            end 
        end
        s::Float64 = 0
        @inbounds for i in eachindex(counts)
            s += bin_ll(bin_centers[i], bin_widths[i], counts[i])
        end
        return s
    end 
    optim_result = optimize( log_likelihood, fit.initial_parameters )
end


"""
     llhfit!(fit::FitFunction, h::Histogram)::Nothing

Performs a log-likelihood fit of the model function `fit.model` and the initial parameters `fit.initial_parameters`
on the histogram `h` in the range `fit.fitrange`. The determined parameters are stored in `fit.parameters`.

The likelihood for each individual bin is the Poission distribution.

There are no uncertainty estimations for this fit yet.
"""
function llhfit!(fit::FitFunction, h::Histogram)::Nothing
    optim_result = llhfit(fit, h)
    fit.parameters = optim_result.minimizer
    nothing
end
