export BATHistogramFitFunction
struct BATHistogramFitFunction{T, ND, NP} <: AbstractFitFunction{T, ND, NP}
    model::Function
    fitranges::NTuple{ND, AbstractVector{T}}
    prior
    backend_result::Vector{Any}
    parameters::Vector{Any}

    function BATHistogramFitFunction{T}(model::Function, ndims::Int, nparams::Int) where {T <: AbstractFloat}
        fitranges::NTuple{ndims, Vector{T}} = NTuple{ndims, Vector{T}}( [-Inf, Inf] for idim in 1:ndims)
        return new{T, ndims, nparams}(model, fitranges, missing, [missing], [missing])
    end
    function BATHistogramFitFunction{T}(model::Function, ndims::Int, prior) where {T <: AbstractFloat}
        fitranges::NTuple{ndims, Vector{T}} = NTuple{ndims, Vector{T}}( [-Inf, Inf] for idim in 1:ndims)
        return new{T, ndims, length(prior)}(model, fitranges, prior, [missing], [missing])
    end
end


@recipe function f(ff::BATHistogramFitFunction; npoints = 501, bin_width = 1.0)
    x = collect(range(ff.fitranges[1][1], stop=ff.fitranges[1][2], length=npoints))
    y = bin_width .* (ff.model(x, collect(parameters)))
    linecolor --> :red
    label --> "Fit model with fitted parameters"
    x,y
end


struct HistogramModelLikelihood{H<:Histogram{<:Real, 1}, F<:BATHistogramFitFunction, T <: Real} <: BAT.AbstractDensity
    h::H
    f::F
    weights::Vector{T}
    midpoints::Vector{T}
    bin_volumes::Vector{T}

    function HistogramModelLikelihood(h::H, f::BATHistogramFitFunction) where {H <: Histogram{<:Real, 1}}
        T = get_pricision_type(f)
        new{H, typeof(f), T}( h, f, h.weights,
                        StatsBase.midpoints(h.edges[1]),
                        T[StatsBase.binvolume(h, i) for i in eachindex(h.weights)] )
    end
end

BAT.nparams(mll::HistogramModelLikelihood) = get_nparams(mll.f)

function BAT.density_logval( l::HistogramModelLikelihood, pars)
    log_likelihood::Float64 = 0
    @inbounds for i in eachindex(l.weights)
        expected_counts::Float64 = l.bin_volumes[i] * l.f.model(l.midpoints[i], pars)
        log_likelihood += logpdf(Poisson(expected_counts), l.weights[i])
    end
    if isnan(log_likelihood)

    end
    return log_likelihood
end

struct BATResult
    result
end

function get_fit_backend_result(r::BATResult)
    return r.result
end


function set_parameters!(ff::BATHistogramFitFunction, pars) 
    ff.parameters[1] = pars
end
function _set_fitted_parameters!(ff::BATHistogramFitFunction)
    set_parameters!(ff, get_fit_backend_result(ff)[3].mode)
end

export batfit!
function batfit!(f::BATHistogramFitFunction, h::Histogram{<:Real, 1};
                nsamples::Int = 1000, nchains::Int = 4, pretunesamples = 4000, max_ncycles = 10, BGConvergenceThreshold = 1.1)
    first_bin::Int = !isinf(first(f.fitranges[1])) ? StatsBase.binindex(h, first(f.fitranges[1])) : 1
    last_bin::Int  = !isinf(last( f.fitranges[1])) ? StatsBase.binindex(h, last( f.fitranges[1])) : length(h.weights)
    if first_bin < 1 first_bin = 1 end
    if (last_bin > length(h.weights)) last_bin = length(h.weights) end
    h_sub = Histogram(h.edges[1][first_bin:last_bin], h.weights[first_bin:last_bin-1])
    h_bat = HistogramModelLikelihood(h_sub, f)

    posterior = PosteriorDensity(h_bat, f.prior)
    algorithm = MetropolisHastings()

    tuner_config = ProposalCovTunerConfig(
        λ = 0.5,
        α = 0.15..0.35,
        β = 1.5,
        c = 1e-4..1e2
    )

    convergence_test = BGConvergence(BGConvergenceThreshold)

    init_strategy = MCMCInitStrategy(
        ninit_tries_per_chain = 8..128,
        max_nsamples_pretune = pretunesamples,
        max_nsteps_pretune = pretunesamples,
        max_time_pretune = Inf
    )

    burnin_strategy = MCMCBurninStrategy(
        max_nsamples_per_cycle = pretunesamples,
        max_nsteps_per_cycle = pretunesamples,
        max_time_per_cycle = Inf,
        max_ncycles = max_ncycles
    )

    bat_result = BATResult(rand(    MCMCSpec(algorithm, posterior), 
                                    nsamples, 
                                    nchains,
                                    tuner_config = tuner_config,
                                    convergence_test = convergence_test,
                                    init_strategy = init_strategy,
                                    burnin_strategy = burnin_strategy,
                                    max_nsteps = nsamples,
                                    max_time = Inf,
                                    granularity = 1 ))

    f.backend_result[1] = bat_result
    _set_fitted_parameters!(f)

    return nothing
end

function get_samples_inds_by_chain_id(sampleids::MCMCSampleIDVector, ichain::Int)
    chainids = Int[]
    for s in sampleids
        if !(s.chainid in chainids) push!(chainids, s.chainid) end
    end
    chainid = chainids[ichain]
    inds = findall( sid -> sid.chainid == chainid, sampleids )
    return inds
end

function get_histogram_pdf(samples::DensitySampleVector; nbins = 10, sample_range::AbstractVector{Int} = Int[]) where {N}
    if isempty(sample_range) sample_range = 1:length(samples) end
    pdf = fit(Histogram, Tuple([flatview(samples.params)[ipar, sample_range] for ipar in 1:length(samples[1])]), 
                FrequencyWeights(samples.weight[sample_range]), closed=:left, nbins=nbins)
    return normalize(pdf)
end

function get_marginalized_pdf(samples::DensitySampleVector, pars::NTuple{N, Int}; nbins = 100, sample_range::AbstractVector{Int} = Int[]) where {N}
    if isempty(sample_range) sample_range = 1:length(samples) end
    pdf = fit(Histogram, Tuple([flatview(samples.params)[ipar, sample_range] for ipar in pars]), 
                FrequencyWeights(samples.weight[sample_range]), closed=:left, nbins=nbins)
    return normalize(pdf)
end
get_marginalized_pdf(samples::DensitySampleVector, ipar::Int; nbins = 100, sample_range::AbstractVector{Int} = Int[]) =
    get_marginalized_pdf(samples, (ipar,); nbins = nbins, sample_range = sample_range)

get_marginalized_pdf(samples::DensitySampleVector, xpar::Int, ypar::Int; nbins = 100, sample_range::AbstractVector{Int} = Int[]) =
    get_marginalized_pdf(samples, (xpar, ypar); nbins = nbins, sample_range = sample_range)

get_marginalized_pdf(samples::DensitySampleVector, var::ShapesOfVariables.VariableDataAccessor; nbins = 100, sample_range::AbstractVector{Int} = Int[]) =
    get_marginalized_pdf(samples, Tuple([var.offset + i for i in 1:var.len ]); nbins = nbins, sample_range = sample_range)

function calculate_localmode(h::StatsBase.Histogram{<:Real, N}) where {N}
    cartesian_inds = CartesianIndices(h.weights)[findmax(h.weights)[2]]
    return [StatsBase.midpoints(h.edges[i])[cartesian_inds[i]] for i in 1:N]
end

function get_marginalized_pdf(f::BATHistogramFitFunction, ipar::Int; nbins = 100)
    @assert typeof(f.backend_result[1]) == BATResult "BATHistogramFitFunction `$f` does not have a `BATResult`."
    samples = f.backend_result[1].result[1]
    pdf = fit(Histogram, flatview(samples.params)[ipar, :], FrequencyWeights(samples.weight), nbins = nbins, closed = :left)
    return normalize(pdf)
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

