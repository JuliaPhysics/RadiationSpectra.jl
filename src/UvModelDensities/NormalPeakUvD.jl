
struct NormalPeakUvD{T} <: UvModelDensity{T}
    # Signal (Normal distribution with area A)
    A::T # Peak area
    UvNormal::Normal{T} # Normal distribution of the peak
    # background (Polynomial of order N)
    bg_poly::Polynomial{2,T}
    bgxl::T # x-value of the left background level
    bgxr::T # x-value of the right background level
    bgll::T # left background level 
    bglr::T # right background level 
end

function evaluate(d::NormalPeakUvD{T}, x::Real) where {T}
    return d.A * pdf(d.UvNormal, x) + d.bg_poly(x)
end

function NormalPeakUvD(nt::NamedTuple{(:A,:μ,:σ,:bgxl,:bgxr,:bgll,:bglr)})
    T = promote_type(typeof.(values(nt))...)
    poly_coeffs = poly_coeffs_from_levels(T(nt.bgxl), T(nt.bgxr), T(nt.bgll), T(nt.bglr))
    NormalPeakUvD(T(nt.A), Normal(T(nt.μ), T(nt.σ)), Polynomial(poly_coeffs), T(nt.bgxl), T(nt.bgxr), T(nt.bgll), T(nt.bglr))
end

function poly_coeffs_from_levels(bgxl, bgxr, bgll, bglr) 
    slope = (bglr - bgll) / (bgxr - bgxl) # linear slope
    mid_offset = bgll + slope * 0.5*(bgxr - bgxl) # offset
    offset = mid_offset - (bgxr + bgxl)/2 * slope
    offset, slope
end

function initial_parameter_guess(d::Type{NormalPeakUvD}, h::StatsBase.Histogram{<:Any, 1}, ::Type{T} = Float64) where {T}
    nh = normalize(h)
    # estimate background parameters through side windows
    bgxl = T(h.edges[1][1])
    bgxr = T(h.edges[1][end])
    Δx_side = (bgxr - bgxl) / 8
    bgll = zero(T)
    bglr = zero(T)
    widths = diff(h.edges[1])
    inv_widths = inv.(widths)
    bgwl = zero(T)
    bgwr = zero(T)
    for i in 1:StatsBase.binindex(h, bgxl + Δx_side)
        bgll += h.weights[i] * inv_widths[i]
        bgwl += widths[i]
    end
    for i in StatsBase.binindex(h, bgxr - Δx_side):length(h.weights)
        bglr += h.weights[i] * inv_widths[i]
        bgwr += widths[i]
    end
    # estimate peak parameters
    μ::T = EmpiricalDistributions._mean(nh)[1]
    σ::T = sqrt(EmpiricalDistributions._var(nh, (μ,))[1])
    A::T = sum(h.weights) - (bgll*bgwl + bglr*bglr)/2
    if A < 0 A = sum(h.weights) end

    parshape = NamedTupleShape(
        A = ScalarShape{T}(),
        μ = ScalarShape{T}(),
        σ = ScalarShape{T}(),
        bgxl = ConstValueShape(bgxl),
        bgxr = ConstValueShape(bgxr),
        bgll = ScalarShape{T}(),
        bglr = ScalarShape{T}()
    )
    p0 = (
        A = A,
        μ = μ,
        σ = σ,
        bgxl = bgxl,
        bgxr = bgxr,
        bgll = bgll,
        bglr = bglr
    )
    lower_bounds = (
        A = T(0),
        μ = μ - σ,
        σ = T(σ/10),
        bgxl = T(bgxl - minimum(widths)),
        bgxr = T(bgxr - minimum(widths)),
        bgll = T(0),
        bglr = T(0)
    )
    upper_bounds = (
        A = T(2*sum(h.weights)),
        μ = μ + σ,
        σ = 3σ,
        bgxl = T(bgxl + minimum(widths)),
        bgxr = T(bgxr + minimum(widths)),
        bgll = T(2bgll),
        bglr = T(2bglr)
    )
    bounds = (
        A = lower_bounds.A..upper_bounds.A,
        μ = lower_bounds.μ..upper_bounds.μ,
        σ = lower_bounds.σ..upper_bounds.σ,
        bgxl = ConstValueDist(bgxl),
        bgxr = ConstValueDist(bgxr),
        bgll = lower_bounds.bgll..upper_bounds.bgll,
        bglr = lower_bounds.bglr..upper_bounds.bglr
    )
    p0, lower_bounds, upper_bounds, parshape, bounds
end

opt_fit(::Type{NormalPeakUvD}, h::Histogram{<:Any, 1}, ::Type{T} = Float64) where {T} =
    opt_fit(RadiationSpectra.NormalPeakUvD, h, initial_parameter_guess(RadiationSpectra.NormalPeakUvD, h, T)[1:4]...) 


 
