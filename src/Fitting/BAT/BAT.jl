export BATResult
export batfit!
export get_marginalized_pdf
export get_histogram_pdf
export get_samples_inds_by_chain_id

struct HistogramModelLikelihood{H<:Histogram{<:Real, 1}, F<:FitFunction, T <: Real} <: BAT.AbstractDensity
    h::H
    f::F
    weights::Vector{T}
    midpoints::Vector{T}
    bin_volumes::Vector{T}
    logabsgamma::Vector{T}
end

function HistogramModelLikelihood(h::H, f::FitFunction) where {H <: Histogram{<:Real, 1}}
    T = get_pricision_type(f)
    HistogramModelLikelihood{H, typeof(f), T}( h, f, h.weights,
                    StatsBase.midpoints(h.edges[1]),
                    T[StatsBase.binvolume(h, i) for i in eachindex(h.weights)],
                    T[logabsgamma(w + 1)[1] for w in h.weights]  )
end

function log_pdf_poisson(λ::T, k::U, logabsgamma::T) where {T<:Real,U<:Real}
    R = float(promote_type(T,U))
    if λ >= 0 && k >= 0 && isinteger(k)
        result = (iszero(k) ? R(k) : R(log(λ) * k)) - λ - logabsgamma
        R(result)
    else
        R(-Inf)
    end
end

function BAT.eval_logval_unchecked(l::HistogramModelLikelihood{H, F, T}, pars) where {H, F, T}
    log_likelihood::Float64 = 0
    @inbounds for i in eachindex(l.weights)
        expected_counts::Float64 = l.bin_volumes[i] * l.f.model(l.midpoints[i], pars)
        if isnan(expected_counts) || expected_counts < 0 
            expected_counts = T(Inf) # This should be removed. The model function, `l.f.model`, should only return values >= 0 
        end
        log_likelihood += log_pdf_poisson(expected_counts, l.weights[i], l.logabsgamma[i])
    end
    return log_likelihood
end


mutable struct BATResult{B, P} 
    bat_result::B
    prior::P

    BATResult(b::B) where B = new{B, Missing}(b, missing)
    BATResult(b::B, p::P) where {B, P} = new{B, P}(b, p)
end


function batfit!(f::FitFunction{T, ND, NP}, h::Histogram{<:Real, 1};
                rng_seed = Philox4x((123, 456)),
                prior = missing,
                nsamples::Int = 10^5, 
                nchains::Int = 4, 
                pretunesamples::Int = length(f.parameter_bounds) * 1000, 
                max_ncycles::Int = 30, 
                tuning_r = 0.5,
                BGConvergenceThreshold::Real = sqrt(length(f.parameter_bounds))) where {T, ND, NP}
    first_bin = !isinf(first(f.fitranges[1])) ? StatsBase.binindex(h, first(f.fitranges[1])) : 1
    last_bin  = !isinf(last( f.fitranges[1])) ? StatsBase.binindex(h, last( f.fitranges[1])) : length(h.weights)
    if first_bin < 1 first_bin = 1 end
    if (last_bin > length(h.weights)) last_bin = length(h.weights) end
    h_sub = Histogram(h.edges[1][first_bin:last_bin], h.weights[first_bin:last_bin-1])
    h_bat = HistogramModelLikelihood(h_sub, f)

    if ismissing(prior)
        prior_dists = if !ismissing(f.parameter_bounds)
            Tuple(map(b -> Uniform(b.left, b.right), f.parameter_bounds))
        else
            error("`batfit!` needs a `FitFunction` with parameter_bounds.")
        end

        prior = BAT.NamedTupleDist( (;zip( f.parameter_names, prior_dists)...))
        prior = BAT.DistributionDensity(prior)#, bounds_type = BAT.reflective_bounds)
    end

    posterior = BAT.PosteriorDensity(h_bat, prior)

    tuning = BAT.AdaptiveMHTuning(
        r = tuning_r,
        λ = 0.5,
        α = 0.15..0.35,
        β = 1.5,
        c = 1e-4..1e2
    )

    mcmcalgo = BAT.MetropolisHastings(
        weighting = BAT.RepetitionWeighting(),
        tuning = tuning
    )

    convergence = BAT.BrooksGelmanConvergence(
        threshold = T(BGConvergenceThreshold),
        corrected = false
    )

    init = BAT.MCMCChainPoolInit(
        init_tries_per_chain = 8..128,
        nsteps_init = pretunesamples * 10,
    )

    burnin = BAT.MCMCMultiCycleBurnin( 
        nsteps_per_cycle = pretunesamples * 10, 
        max_ncycles = max_ncycles
    )

    bat_result = BATResult(
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
        ),
        prior
    )

    f.backend_result = (bat_result, )
    try 
        _set_fitted_parameters!(f, mode(RadiationSpectra.get_bat_result(f))[])
    catch err
        @warn "Determination of mode failed"
    end
    f
end

### Utilties

function rs_write(fn::AbstractString, b::BATResult)
    BAT.bat_write(fn, unshape(b))   
    ser_bfn = begin
        bn = basename(fn)
        if occursin(".", bn)
            bn[1:first(findlast(".", bn))-1] * ".ser"
        else
            bn * ".ser"
        end
    end
    ser_fn = joinpath(dirname(fn), ser_bfn)
    serialize(ser_fn, b.prior)
    nothing
end
function rs_read(fn::AbstractString)
    ser_bfn = begin
        bn = basename(fn)
        if occursin(".", bn)
            bn[1:first(findlast(".", bn))-1] * ".ser"
        else
            bn * ".ser"
        end
    end
    ser_fn = joinpath(dirname(fn), ser_bfn)
    if !isfile(ser_fn) error("Prior was not serialized to file.") end
    prior = deserialize(ser_fn)
    s = varshape(prior)
    bat_samples = s.(BAT.bat_read(fn).result)
    return BATResult((result = bat_samples,), prior)
end


get_fit_backend_result(fr::Tuple) = fr

get_bat_result(ff::FitFunction)::BAT.DensitySampleVector = 
    get_fit_backend_result(ff.backend_result)[1].bat_result.result


function BAT.NamedTupleShape(f::FitFunction)
    ts = Tuple(map(b -> ScalarShape{Real}(), f.parameter_names))
    return BAT.NamedTupleShape(NamedTuple{Tuple(f.parameter_names)}(ts))
end

function _set_fitted_parameters!(ff::FitFunction{T}, fitted_parameters::ShapedAsNT)::Nothing where {T <: AbstractFloat}
    _set_fitted_parameters!(ff, vcat(values(fitted_parameters[])...))
end


get_samples_inds_by_chain_id(ff::FitFunction, ichain::Int = 1) = 
    get_samples_inds_by_chain_id(get_bat_result(ff), ichain)

function get_samples_inds_by_chain_id(bat_result, ichain::Int)
    chainids = Int[]
    sampleids = bat_result.info
    for s in sampleids
        if !(s.chainid in chainids) push!(chainids, s.chainid) end
    end
    chainid = chainids[ichain]
    inds = findall( sid -> sid.chainid == chainid, sampleids )
    return inds
end

function _get_standard_deviations(f::FitFunction{T}, fr::Tuple) where {T}
    uncertainties = zeros(T, length(f.fitted_parameters))
    for i in eachindex(f.fitted_parameters)  
        h = get_marginalized_pdf(BAT.unshaped.(get_bat_result(f)), i, nbins = 100)
        d = fit(Normal, h.edges[1], h.weights)
        uncertainties[i] = d.σ
    end
    return uncertainties
end

unshape(samples::BAT.DensitySampleVector) =
    samples.v isa ShapedAsNTArray ? ValueShapes.unshaped.(samples) : samples
unshape(b::BATResult) =
    unshape(b.bat_result[1])


function get_marginalized_pdf(samples::BAT.DensitySampleVector, pars::NTuple{N, Int}; nbins = 100, sample_range::AbstractVector{Int} = Int[]) where {N}
    if isempty(sample_range) sample_range = 1:length(samples) end
    unshaped_samples = unshape(samples)
    pdf = fit(Histogram, Tuple([flatview(unshaped_samples.v)[ipar, sample_range] for ipar in pars]), 
                FrequencyWeights(unshaped_samples.weight[sample_range]), closed=:left, nbins=nbins)
    return normalize(pdf)
end

get_marginalized_pdf(br::BATResult, args...; kwargs...) =
    get_marginalized_pdf(br.bat_result.result, args...; kwargs...)

get_marginalized_pdf(f::FitFunction, args...; nbins = 100, sample_range::AbstractVector{Int} = Int[]) =
    get_marginalized_pdf(get_bat_result(f), args..., nbins = nbins, sample_range = sample_range)

get_marginalized_pdf(samples::BAT.DensitySampleVector, ipar::Int; nbins = 100, sample_range::AbstractVector{Int} = Int[]) =
    get_marginalized_pdf(samples, (ipar,); nbins = nbins, sample_range = sample_range)
    
get_marginalized_pdf(samples::BAT.DensitySampleVector, var::ValueShapes.ValueAccessor; nbins = 100, sample_range::AbstractVector{Int} = Int[]) =
    get_marginalized_pdf(samples, Tuple([var.offset + i for i in 1:var.len ]); nbins = nbins, sample_range = sample_range)
    

function get_histogram_pdf(samples::BAT.DensitySampleVector; nbins = 10, sample_range::AbstractVector{Int} = Int[]) where {N}
    unshaped_samples = unshape(samples)
    n_params::Int = size(flatview(unshaped_samples.v), 1)
    get_marginalized_pdf(unshaped_samples, Tuple(1:n_params), nbins = nbins, sample_range = sample_range)
end
get_histogram_pdf(f::FitFunction, args...; kwargs...) = 
    get_histogram_pdf(get_bat_result(f), args...; kwargs...)



function calculate_localmode(h::StatsBase.Histogram{<:Real, N}) where {N}
    cartesian_inds = CartesianIndices(h.weights)[findmax(h.weights)[2]]
    return [StatsBase.midpoints(h.edges[i])[cartesian_inds[i]] for i in 1:N]
end


function central_intervals(marginalized_pdf::Histogram{<:Real, 1}, intervals = BAT.default_credibilities)
    sort!(intervals)
    splitted_histograms = reverse(BAT.split_central(marginalized_pdf, intervals)[1])
    r = Interval[]
    for ipar in 1:length(intervals)
        ifirst = findfirst(w -> w > 0, splitted_histograms[ipar].weights)
        ilast  = findlast( w -> w > 0, splitted_histograms[ipar].weights)
        int = splitted_histograms[ipar].edges[1][ifirst]..splitted_histograms[ipar].edges[1][ilast+1]
        push!(r, int)
    end
    return r
end
function smallest_intervals(marginalized_pdf::Histogram{<:Real, 1}, intervals = BAT.default_credibilities)
    sort!(intervals)
    splitted_histograms = reverse(BAT.split_smallest(marginalized_pdf, intervals)[1])
    r = Interval[]
    for ipar in 1:length(intervals)
        ifirst = findfirst(w -> w > 0, splitted_histograms[ipar].weights)
        ilast  = findlast( w -> w > 0, splitted_histograms[ipar].weights)
        int = splitted_histograms[ipar].edges[1][ifirst]..splitted_histograms[ipar].edges[1][ilast+1]
        push!(r, int)
    end
    return r
end

@userplot Plot_Marginalized_PDF
@recipe function f(h::Plot_Marginalized_PDF; interval_type = :central)
    if (length(h.args) < 1) || !(typeof(h.args[1]) <: Histogram{<:Real, 1}) 
        error("Wrong arguments. Got: $(typeof(h.args))")
    end
    hmpdf = h.args[1]
    intervals = BAT.default_credibilities

    @series begin
        seriestype := :step
        label --> "Marginalized PDF"
        hmpdf
    end

    h_ints = if interval_type == :smallest
        reverse(BAT.split_smallest(hmpdf, BAT.default_credibilities)[1])
    elseif interval_type == :central
        reverse(BAT.split_central(hmpdf, BAT.default_credibilities)[1])
    else
        error("`interval_type` must be either `:smallest` or `:central`")
    end
    for i in length(intervals):-1:1
        @series begin
            label := """$((interval_type == :smallest ? "SI" : "CI")): $(BAT.default_credibilities[i])"""
            fillcolor := BAT.default_colors[i] 
            h_ints[i]
        end
    end
end

