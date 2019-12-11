struct HistogramModelLikelihood{H<:Histogram{<:Real, 1}, F<:FitFunction, T <: Real} <: BAT.AbstractDensity
    h::H
    f::F
    weights::Vector{T}
    midpoints::Vector{T}
    bin_volumes::Vector{T}

    function HistogramModelLikelihood(h::H, f::FitFunction) where {H <: Histogram{<:Real, 1}}
        T = get_pricision_type(f)
        new{H, typeof(f), T}( h, f, h.weights,
                        StatsBase.midpoints(h.edges[1]),
                        T[StatsBase.binvolume(h, i) for i in eachindex(h.weights)] )
    end
end

# BAT.nparams(mll::HistogramModelLikelihood) = get_nparams(mll.f)

function BAT.density_logval( l::HistogramModelLikelihood{H, F, T}, pars) where {H, F, T}
    log_likelihood::Float64 = 0
    @inbounds for i in eachindex(l.weights)
        expected_counts::Float64 = l.bin_volumes[i] * l.f.model(l.midpoints[i], pars)
        if isnan(expected_counts) || expected_counts < 0 
            expected_counts = T(Inf)
        end
        log_likelihood += logpdf(Poisson( expected_counts ), l.weights[i])
    end
    return log_likelihood
end


export batfit!
function batfit!(f::FitFunction{T, ND, NP}, h::Histogram{<:Real, 1};
                nsamples::Int = 10^5, 
                nchains::Int = 4, 
                pretunesamples::Int = length(f.parameter_bounds) * 1000, 
                max_ncycles::Int = 30, 
                BGConvergenceThreshold::Real = sqrt(length(f.parameter_bounds))) where {T, ND, NP}
    first_bin = !isinf(first(f.fitranges[1])) ? StatsBase.binindex(h, first(f.fitranges[1])) : 1
    last_bin  = !isinf(last( f.fitranges[1])) ? StatsBase.binindex(h, last( f.fitranges[1])) : length(h.weights)
    if first_bin < 1 first_bin = 1 end
    if (last_bin > length(h.weights)) last_bin = length(h.weights) end
    h_sub = Histogram(h.edges[1][first_bin:last_bin], h.weights[first_bin:last_bin-1])
    h_bat = HistogramModelLikelihood(h_sub, f)

    prior_dists = if !ismissing(f.parameter_bounds)
        Tuple(map(b -> Uniform(b.left, b.right), f.parameter_bounds))
    else
        error("`batfit!` needs a `FitFunction` with parameter_bounds.")
    end

    prior = BAT.NamedTupleDist( (;zip( f.parameter_names, prior_dists)...))

    posterior = BAT.PosteriorDensity(h_bat, prior)
    algorithm = BAT.MetropolisHastings()

    tuning = BAT.AdaptiveMetropolisTuning(
        λ = 0.5,
        α = 0.15..0.35,
        β = 1.5,
        c = 1e-4..1e2
    )

    convergence = BAT.BrooksGelmanConvergence(
        threshold = T(BGConvergenceThreshold),
        corrected = false
    )

    init = BAT.MCMCInitStrategy(
        ninit_tries_per_chain = 8..128,
        max_nsamples_pretune = pretunesamples,
        max_nsteps_pretune = pretunesamples * 10,
        max_time_pretune = Inf
    )

    burnin = BAT.MCMCBurninStrategy(
        max_nsamples_per_cycle = 10^3,
        max_nsteps_per_cycle = 10^4,
        max_time_per_cycle = Inf,
        max_ncycles = max_ncycles
    )

    samples, stats = BAT.bat_sample(
        posterior, (nsamples, nchains), algorithm,
        max_nsteps = 10 * nsamples,
        max_time = Inf,
        tuning = tuning,
        init = init,
        burnin = burnin,
        convergence = convergence,
        strict = false,
        filter = true
    )

    f.backend_result = (samples, stats) 
   _set_fitted_parameters!(f, stats.mode)
    f
end


function get_samples_inds_by_chain_id(samples::BAT.PosteriorSampleVector, ichain::Int)
    chainids = Int[]
    sampleids = samples.info
    for s in sampleids
        if !(s.chainid in chainids) push!(chainids, s.chainid) end
    end
    chainid = chainids[ichain]
    inds = findall( sid -> sid.chainid == chainid, sampleids )
    return inds
end

function get_histogram_pdf(samples::BAT.PosteriorSampleVector; nbins = 10, sample_range::AbstractVector{Int} = Int[]) where {N}
    if isempty(sample_range) sample_range = 1:length(samples) end
    n_params::Int = length(samples[1].params)
    pdf = fit(Histogram, Tuple([flatview(samples.params)[ipar, sample_range] for ipar in 1:n_params]), 
                FrequencyWeights(samples.weight[sample_range]), closed=:left, nbins=nbins)
    return normalize(pdf)
end

function get_marginalized_pdf(samples::BAT.PosteriorSampleVector, pars::NTuple{N, Int}; nbins = 100, sample_range::AbstractVector{Int} = Int[]) where {N}
    if isempty(sample_range) sample_range = 1:length(samples) end
    pdf = fit(Histogram, Tuple([flatview(samples.params)[ipar, sample_range] for ipar in pars]), 
                FrequencyWeights(samples.weight[sample_range]), closed=:left, nbins=nbins)
    return normalize(pdf)
end
get_marginalized_pdf(samples::BAT.PosteriorSampleVector, ipar::Int; nbins = 100, sample_range::AbstractVector{Int} = Int[]) =
    get_marginalized_pdf(samples, (ipar,); nbins = nbins, sample_range = sample_range)

get_marginalized_pdf(samples::BAT.PosteriorSampleVector, xpar::Int, ypar::Int; nbins = 100, sample_range::AbstractVector{Int} = Int[]) =
    get_marginalized_pdf(samples, (xpar, ypar); nbins = nbins, sample_range = sample_range)

get_marginalized_pdf(samples::BAT.PosteriorSampleVector, var::ValueShapes.ValueAccessor; nbins = 100, sample_range::AbstractVector{Int} = Int[]) =
    get_marginalized_pdf(samples, Tuple([var.offset + i for i in 1:var.len ]); nbins = nbins, sample_range = sample_range)

function calculate_localmode(h::StatsBase.Histogram{<:Real, N}) where {N}
    cartesian_inds = CartesianIndices(h.weights)[findmax(h.weights)[2]]
    return [StatsBase.midpoints(h.edges[i])[cartesian_inds[i]] for i in 1:N]
end

function get_marginalized_pdf(f::FitFunction, ipar::Int; nbins = 100)
    return get_marginalized_pdf(f.backend_result[1], ipar, nbins = nbins)
end

function central_intervals(marginalized_pdf::Histogram{<:Real, 1}, intervals = BAT.standard_confidence_vals)
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

@userplot Plot_Marginalized_PDF
@recipe function f(h::Plot_Marginalized_PDF)
    if (length(h.args) < 1) || !(typeof(h.args[1]) <: Histogram{<:Real, 1}) 
        error("Wrong arguments. Got: $(typeof(h.args))")
    end
    hmpdf = h.args[1]
    intervals = BAT.standard_confidence_vals

    @series begin
        seriestype := :step
        label --> "Marginalized PDF"
        hmpdf
    end

    h_ints = reverse(BAT.split_central(hmpdf, BAT.standard_confidence_vals)[1])
    for i in length(intervals):-1:1
        @series begin
            label := "CI: $(BAT.standard_confidence_vals[i])"
            fillcolor := BAT.standard_colors[i] 
            h_ints[i]
        end
    end
end


get_fit_backend_result(fr::Tuple{<:Any,<:Any}) = fr

function _get_standard_deviations(f::FitFunction{T}, fr::Tuple{<:AbstractArray{<:BAT.PosteriorSample}, <:Any}) where {T}
    uncertainties = zeros(T, length(f.fitted_parameters))
    for i in eachindex(f.fitted_parameters)  
        h = get_marginalized_pdf(fr[1], i, nbins = 100)
        d = fit(Normal, h.edges[1], h.weights)
        uncertainties[i] = d.σ
    end
    return uncertainties
end


