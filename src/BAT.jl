abstract type BATFitBackend <: FitBackend end
backend_type(::Val{:BAT}) = BATFitBackend
fit(::Type{BATFitBackend}, args...; kwargs...) = bat_fit(args...; kwargs...)

function bat_fit(UvD::Type{<:UvSpectrumDensity}, h::Histogram{<:Any, 1};
            bounds = initial_parameter_guess(UvD, h)[end],
            rng_seed = Philox4x((321, 456)),
            prior = missing,
            nsamples::Int = 10^5, 
            nchains::Int = 4, 
            nsteps_init::Int = 2000, 
            nsteps_per_cycle::Int = 2000, 
            strict = true,
            mcmcalgo = BAT.HamiltonianMC(),
            convergence = missing, 
            init = missing, 
            burnin = missing, 
            max_ncycles::Int = 1, 
            threshold::Real = 1.1
)
    hp = HistLLHPrecalulations(float(h))
    sd = SpectrumDensity{typeof(hp), UvD}(hp)
    posterior = BAT.PosteriorDensity(sd, ismissing(prior) ? BAT.NamedTupleDist(bounds) : prior )
    posterior_transform = BAT.bat_transform(BAT.PriorToGaussian(), posterior, BAT.PriorSubstitution()).result
    trafo = BAT.trafoof(posterior_transform.likelihood) 

    convergence = if ismissing(convergence)
        BAT.BrooksGelmanConvergence(
            threshold = threshold,
            corrected = false
        )
    else
        convergence
    end

    init = if ismissing(init)
        BAT.MCMCChainPoolInit(
            init_tries_per_chain = 8..128,
            nsteps_init = nsteps_init,
        )
    else
        init
    end

    burnin = if ismissing(burnin)
        BAT.MCMCMultiCycleBurnin( 
            nsteps_per_cycle = nsteps_per_cycle, 
            max_ncycles = max_ncycles
        )
    else
        burnin
    end

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
    return UvD(mode(bat_samples)), bat_samples
end

