abstract type BATFitBackend <: FitBackend end
backend_type(::Val{:BAT}) = BATFitBackend
fit(::Type{BATFitBackend}, args...; kwargs...) = bat_fit(args...; kwargs...)

struct BATHistLLHPrecalulations{H, DT} <: BAT.AbstractDensity
    d::H
end
BATHistLLHPrecalulations(d::H, ::Type{DT}) where {H, DT} = BATHistLLHPrecalulations{H, DT}(d)

function BAT.eval_logval_unchecked(d::BATHistLLHPrecalulations{H, DT}, pars) where {H, DT}
    log_likelihood = zero(eltype(pars))
    model_dist = DT(pars)
    @inbounds for i in eachindex(d.d.weights)
        expected_counts = d.d.volumes[i] * evaluate(model_dist, d.d.midpoints[i])
        log_likelihood += log_pdf_poisson(expected_counts, d.d.weights[i], d.d.logabsgamma[i])
    end
    return log_likelihood
end

function bat_fit(UvD::Type{<:UvModelDensity}, h::Histogram{<:Any, 1};
            bounds = initial_parameter_guess(UvD, h)[end],
            rng_seed = Philox4x((321, 456)),
            prior = missing,
            nsamples::Int = 10^5, 
            nchains::Int = 4, 
            nsteps_init::Int = 2000, 
            nsteps_per_cycle::Int = 2000, 
            strict = true,
            max_ncycles::Int = 1, 
            threshold::Real = 1.1) 
    hp = BATHistLLHPrecalulations(HistLLHPrecalulations(h), UvD)
    posterior = BAT.PosteriorDensity(hp, ismissing(prior) ? BAT.DistributionDensity(BAT.NamedTupleDist(bounds)) : prior )
    posterior_transform = BAT.bat_transform(BAT.PriorToGaussian(), posterior, BAT.PriorSubstitution()).result
    trafo = BAT.trafoof(posterior_transform.likelihood) 

    mcmcalgo = BAT.HamiltonianMC()

    convergence = BAT.BrooksGelmanConvergence(
        threshold = threshold,
        corrected = false
    )

    init = BAT.MCMCChainPoolInit(
        init_tries_per_chain = 8..128,
        nsteps_init = nsteps_init,
    )

    burnin = BAT.MCMCMultiCycleBurnin( 
        nsteps_per_cycle = nsteps_per_cycle, 
        max_ncycles = max_ncycles
    )

    bat_samples = inv(trafo).(
            BAT.bat_sample( 
                rng_seed, posterior_transform,
                BAT.MCMCSampling(
                    mcalg = mcmcalgo,
                    nchains = nchains,
                    nsteps = nsamples,
                    init = init,
                    burnin = burnin,
                    strict = strict,
                    convergence = convergence,
                )
            ).result
        )
    return UvD(mode(bat_samples)[1]), bat_samples
end

