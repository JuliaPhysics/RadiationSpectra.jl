struct BATHistLLHPrecalulations{H, DT} <: BAT.AbstractDensity
    d::H
end
BATHistLLHPrecalulations(d::H, ::Type{DT}) where {H, DT} = BATHistLLHPrecalulations{H, DT}(d)

function BAT.eval_logval_unchecked(d::BATHistLLHPrecalulations{H, DT}, pars) where {H, DT}
    T = eltype(pars)
    log_likelihood::T = zero(T)
    model_dist = DT(pars)
    @inbounds for i in eachindex(d.d.weights)
        expected_counts = d.d.volumes[i] * evaluate(model_dist, d.d.midpoints[i])
        log_likelihood += log_pdf_poisson(expected_counts, d.d.weights[i], d.d.logabsgamma[i])
    end
    return log_likelihood
end

function batfit(DT::Type{<:UvModelDensity}, h::Histogram{<:Any, 1}, prior::BAT.DistributionDensity;
        rng_seed = Philox4x((123, 456)),
        nsteps::Int = 10^4, 
        nchains::Int = 4, 
        pretunesamples::Int = 100, 
        max_ncycles::Int = 30, 
        BGConvergenceThreshold::Real = sqrt(totalndof(prior)))
    hp = BATHistLLHPrecalulations(HistLLHPrecalulations(h), NormalPeakUvD)
    posterior = BAT.PosteriorDensity(hp, prior)
    posterior = BAT.bat_transform(BAT.PriorToGaussian(), posterior, BAT.PriorSubstitution()).result
    trafo = BAT.trafoof(posterior.likelihood) 

    mcmcalgo = BAT.HamiltonianMC()
    convergence = BAT.BrooksGelmanConvergence(
        threshold = BGConvergenceThreshold,
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

    bat_samples = inv(trafo).(
        BAT.bat_sample( 
            rng_seed, posterior,
            BAT.MCMCSampling(
                mcalg = mcmcalgo,
                nchains = nchains,
                nsteps = nsteps,
                init = init,
                burnin = burnin,
                convergence = convergence,
            )
        ).result
    )
    return DT(mode(bat_samples)[]), bat_samples
end

function batfit(DT::Type{<:UvModelDensity}, h::Histogram{<:Any, 1}, bounds::NamedTuple; kwargs...)
    prior = BAT.DistributionDensity(BAT.NamedTupleDist( bounds ) )
    batfit(NormalPeakUvD, h, prior; kwargs...) 
end

function fit(d::Type{<:AbstractModelDensity}, h::Histogram, backend::Symbol = :Optim, args...; kwargs...)
    if backend == :Optim
        fit(d, h, args...; kwargs...)
    elseif backend == :BAT
        batfit(d, h, args...; kwargs...)
    else
        error("backend must be either `:Optim` (Likelihood optimizer fit via Optim.jl) or `:BAT` (Bayesian fit via BAT.jl). It was $backend.")
    end
end



function batfit(::Type{NormalPeakUvD}, h::Histogram{<:Any, 1}; kwargs...)
    bounds = initial_parameter_guess(NormalPeakUvD{Float64}, h)[end]
    batfit(NormalPeakUvD, h, bounds; kwargs...)
end
