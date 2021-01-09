struct BATHistLLHPrecalulations{H, DT} <: BAT.AbstractDensity
    d::H
end
BATHistLLHPrecalulations(d::H, ::Type{DT}) where {H, DT} = BATHistLLHPrecalulations{H, DT}(d)

function BAT.eval_logval_unchecked(d::BATHistLLHPrecalulations{H, DT}, pars) where {H, DT}
    log_likelihood = 0
    model_dist = DT(pars)
    @inbounds for i in eachindex(d.d.weights)
        expected_counts = d.d.volumes[i] * evaluate(model_dist, d.d.midpoints[i])
        log_likelihood += log_pdf_poisson(expected_counts, d.d.weights[i], d.d.logabsgamma[i])
    end
    return log_likelihood
end


function rsbatfit(::Type{NormalPeakUvD}, h::Histogram{<:Any, 1};
            rng_seed = Philox4x((123, 456)),
            prior = missing,
            nsamples::Int = 10^5, 
            nchains::Int = 2, 
            pretunesamples::Int = 20000, 
            max_ncycles::Int = 30, 
            tuning_r = 0.5,
            BGConvergenceThreshold::Real = sqrt(5) )
    T = Float64
    bounds = initial_parameter_guess(h, NormalPeakUvD{T})[end]

    hp = BATHistLLHPrecalulations(HistLLHPrecalulations(h), NormalPeakUvD)
    prior = BAT.DistributionDensity(BAT.NamedTupleDist( bounds ) )
    posterior = BAT.PosteriorDensity(hp, prior)
    posterior = BAT.bat_transform(BAT.PriorToGaussian(), posterior, BAT.PriorSubstitution()).result
    trafo = BAT.trafoof(posterior.likelihood) 

    mcmcalgo = BAT.HamiltonianMC()
    # mcmcalgo = BAT.MetropolisHastings(
    #     weighting = BAT.RepetitionWeighting(),
    #     tuning = tuning
    # )

    convergence = BAT.BrooksGelmanConvergence(
        threshold = T(BGConvergenceThreshold),
        corrected = false
    )

    init = BAT.MCMCChainPoolInit(
        init_tries_per_chain = 8..128,
        nsteps_init = pretunesamples,
    )

    burnin = BAT.MCMCMultiCycleBurnin( 
        nsteps_per_cycle = pretunesamples, 
        max_ncycles = max_ncycles
    )

    # bat_samples = inv(trafo).(BAT.bat_sample( 
    bat_samples = inv(trafo).(
            BAT.bat_sample( 
                rng_seed, posterior,
                BAT.MCMCSampling(
                    mcalg = mcmcalgo,
                    nchains = nchains,
                    nsteps = nsamples,
                    init = init,
                    burnin = burnin,
                    convergence = convergence,
                )
            ).result
        )
    return bat_samples
end

