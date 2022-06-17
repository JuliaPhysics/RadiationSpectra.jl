struct HistLLHPrecalulations{N, WT, MPT, BVT, LGT} 
    weights::Array{WT, N}
    midpoints::MPT
    volumes::BVT
    logabsgamma::Array{LGT, N}
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

function HistLLHPrecalulations(h::Histogram{<:Any, 1})
    T = promote_type(eltype(h.weights), eltype(h.edges[1]))
    @assert all(isinteger, h.weights) """
        Histogram has non-integer weights. 
        This can not be as the Poisson distribution, the underlying likehood for the bins,
        needs discrete numbers.
    """
    return HistLLHPrecalulations( 
        h.weights,
        StatsBase.midpoints(h.edges[1]),
        T[StatsBase.binvolume(h, i) for i in eachindex(h.weights)],
        T[logabsgamma(w + 1)[1] for w in h.weights]
    )
end

struct SpectrumDensity{HP <: HistLLHPrecalulations, D <: UvSpectrumDensity}
    hp::HP
end

@inline DensityInterface.DensityKind(::SpectrumDensity) = IsDensity()
function DensityInterface.logdensityof(object::SpectrumDensity{HP, D}, x) where {HP, D}
    llh = zero(eltype(x))
    for i in eachindex(object.hp.weights)
        llh += log_pdf_poisson(object.hp.volumes[i] * evaluate(D(x), object.hp.midpoints[i]), object.hp.weights[i], object.hp.logabsgamma[i]) 
    end
    llh
end



"""
    # opt_fit(DT::Type{<:UvSpectrumDensity}, h::Histogram{<:Any, 1},
             p0::AbstractVector, lower_bounds::AbstractVector, upper_bounds::AbstractVector)

Maximum Likelihood Estimation Fit of model density `d` on the histogram `h`.     
"""
function opt_fit(DT::Type{<:UvSpectrumDensity}, h::Histogram{<:Any, 1},  
             p0::AbstractVector, lower_bounds::AbstractVector, upper_bounds::AbstractVector,
             parshape::ValueShapes.AbstractValueShape = valshape(p0)) 
    hp = HistLLHPrecalulations(h)
    function f(p::AbstractVector{T}, hp=hp, ps=parshape) where {T} 
        x = _par_in_input_form(ps(p))
        d = SpectrumDensity{typeof(hp), DT}(hp)
        return -DensityInterface.logdensityof(d, x)
    end
    opt_result = Optim.optimize( f, promote(lower_bounds, upper_bounds, p0)..., 
        Fminbox(BFGS()); autodiff=:forward )
    # return opt_result.minimizer, opt_result
    return DT(_par_in_input_form(parshape(opt_result.minimizer))), opt_result
end


"""
    # opt_fit(DT::Type{<:UvSpectrumDensity}, h::Histogram{<:Any, 1},
             p0::NamedTuple, lower_bounds::NamedTuple, upper_bounds::NamedTuple)

Maximum Likelihood Estimation Fit of model density `d` on the histogram `h`.     
"""
function opt_fit(DT::Type{<:UvSpectrumDensity}, h::Histogram{<:Any, 1},
             p0::NamedTuple, 
             lower_bounds::NamedTuple, 
             upper_bounds::NamedTuple, 
             parshape::NamedTupleShape = valshape(p0))
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



