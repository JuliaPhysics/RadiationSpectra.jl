struct HistogramModelLikelihood{H<:Histogram{<:Real, 1}, F<:FitFunction, T <: Real} <: BAT.AbstractDensity
    h::H
    f::F
    weights::Vector{T}
    midpoints::Vector{T}
    bin_volumes::Vector{T}

    function HistogramModelLikelihood(h::H, f::FitFunction) where {H <: Histogram{<:Real, 1}}
        T = get_pricision_type(f)
        new{H, FitFunction, T}( h, f, h.weights,
                        StatsBase.midpoints(h.edges[1]),
                        T[StatsBase.binvolume(h, i) for i in eachindex(h.weights)] )
    end
end
BAT.nparams(mll::HistogramModelLikelihood) = get_nparams(mll.f)
function BAT.unsafe_density_logval( l::HistogramModelLikelihood, pars,
                                    exec_contect::ExecContext)
    log_likelihood::Float64 = 0
    @inbounds for i in eachindex(l.weights)
        expected_counts::Float64 = l.bin_volumes[i] * l.f.model(l.midpoints[i], pars)
        log_likelihood += logpdf(Poisson(expected_counts), l.weights[i])
    end
    return log_likelihood
end

struct BATResult
    result
end

export batfit!
function batfit!(fit::FitFunction, h::Histogram{<:Real, 1};
                nsamples::Int = 1000, nchains::Int = 4)
    h_bat = HistogramModelLikelihood(h, fit)
    par_bounds = HyperRectBounds( fit.parameter_bounds, reflective_bounds)
    prior = ConstDensity( par_bounds, normalize)
    bayesian_model = BayesianModel(h_bat, prior)
    algorithm = MetropolisHastings()

    bat_result = BATResult(rand( MCMCSpec(algorithm, bayesian_model), nsamples, nchains))
    par_mode = bat_result.result[3].mode
    par_mean = bat_result.result[3].param_stats.mean
    par_covar = bat_result.result[3].param_stats.cov

    _set_fitted_parameters!(fit, par_mode)

    fit.backend_result[1] = bat_result

    return nothing
end
