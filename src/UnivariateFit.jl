

function expectation_values(h::StatsBase.Histogram{<:Any, 1}, DT::Type{<:UvModelDensity}, pars) 
    mps = StatsBase.midpoints(h.edges[1])
    volumes = [StatsBase.binvolume(h, i) for i in eachindex(mps)]
    d = DT(pars)
    [volumes[i] * evaluate(d, mps[i]) for i in eachindex(mps)]
end

function log_pdf_poisson(λ::T, k::U, logabsgamma::LAT) where {T, U, LAT}
    R = float(promote_type(T,U))
    if λ >= 0 && k >= 0 && isinteger(k)
        result = (iszero(k) ? R(k) : R(log(λ) * k)) - λ - logabsgamma
        R(result)
    else
        R(-Inf)
    end
end

struct HistLLHPrecalulations{N, WT, MPT, BVT, LGT} 
    weights::Array{WT, N}
    midpoints::MPT
    volumes::BVT
    logabsgamma::Array{LGT, N}
end

function HistLLHPrecalulations(h::Histogram{<:Any, 1})
    T = promote_type(eltype(h.weights), eltype(h.edges[1]))
    HistLLHPrecalulations( 
        h.weights,
        StatsBase.midpoints(h.edges[1]),
        T[StatsBase.binvolume(h, i) for i in eachindex(h.weights)],
        T[logabsgamma(w + 1)[1] for w in h.weights]  )
end

function loglikelihood(h::HistLLHPrecalulations{1}, DT::Type{<:UvModelDensity}, pars::Vector{PT}, parshape)::PT where PT
    log_likelihood = zero(PT)
    d = DT(_par_in_input_form(parshape(pars)))
    @inbounds for i in eachindex(h.weights)
        expected_counts = h.volumes[i] * evaluate(d, h.midpoints[i])
        log_likelihood += log_pdf_poisson(expected_counts, h.weights[i], h.logabsgamma[i])
    end
    return log_likelihood
end
 
"""
    # opt_fit(DT::Type{<:UvModelDensity}, h::Histogram{<:Any, 1},
             p0::AbstractVector, lower_bounds::AbstractVector, upper_bounds::AbstractVector)

Maximum Likelihood Estimation Fit of model density `d` on the histogram `h`.     
"""
function opt_fit(DT::Type{<:UvModelDensity}, h::Histogram{<:Any, 1},  
             p0::AbstractVector, lower_bounds::AbstractVector, upper_bounds::AbstractVector,
             parshape::ValueShapes.AbstractValueShape = valshape(p0)) 
    hp = HistLLHPrecalulations(h)
    f(p::AbstractVector{T}, hp=hp, DT=DT, ps=parshape) where {T} = -loglikelihood(hp, DT, p, ps)::T 
    opt_result = Optim.optimize( f, promote(lower_bounds, upper_bounds, p0)..., 
        Fminbox(BFGS()); autodiff=:forward )
    # return opt_result.minimizer, opt_result
    return DT(_par_in_input_form(parshape(opt_result.minimizer))), opt_result
end


"""
    # opt_fit(DT::Type{<:UvModelDensity}, h::Histogram{<:Any, 1},
             p0::NamedTuple, lower_bounds::NamedTuple, upper_bounds::NamedTuple)

Maximum Likelihood Estimation Fit of model density `d` on the histogram `h`.     
"""
function opt_fit(DT::Type{<:UvModelDensity}, h::Histogram{<:Any, 1},
             p0::NamedTuple, 
             lower_bounds::NamedTuple, 
             upper_bounds::NamedTuple, 
             parshape::AbstractValueShape = valshape(p0))
    # flatten parameter and bounds
    flat_p0 = ValueShapes.unshaped(p0, parshape)
    T = eltype(flat_p0)
    rnt_lower = _remove_const_shapes_from_nt(lower_bounds, parshape)
    rnt_upper = _remove_const_shapes_from_nt(upper_bounds, parshape)
    flat_lower::Vector{T} = ValueShapes.unshaped(rnt_lower, valshape(rnt_lower))
    flat_upper::Vector{T} = ValueShapes.unshaped(rnt_upper, valshape(rnt_upper))
    
    # return parshape(minimizer)[], opt_result
    return opt_fit(DT, Histogram(h.edges[1], T.(h.weights)), flat_p0, flat_lower, flat_upper, parshape)
end



