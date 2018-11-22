"""
    gauss(x::T, p::Vector{T})::T

A Gaussian with 3 parameters:
- `p[1]`: Scale/Amplitude
- `p[2]`: σ
- `p[3]`: μ
"""
function gauss(x::T, p::Vector{T})::T where {T}
    scale::T = p[1]
    σ::T     = p[2]
    μ::T     = p[3]
    return scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2)) 
end
"""
    gauss(x::Vector{T}, p::Vector{T})::Vector{T}

Maps `x` to `gauss(x::T, p::Vector{T})::T`
"""
function gauss(x::Vector{T}, p::Vector{T})::Vector{T} where {T}
    return T[ gauss(v, p) for v in x ]
end

function gauss_plus_first_order_polynom(x::T, p::Vector{T})::T where {T}
    scale::T = p[1]
    σ::T     = p[2]
    μ::T     = p[3]
    cp0::T   = p[4] 
    cp1::T   = p[5]
    return @fastmath scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2)) + cp0 + cp1 * (x - μ)
end

function gauss_plus_first_order_polynom(x::Vector{T}, p::Vector{T})::Vector{T} where {T}
    return T[ gauss_plus_first_order_polynom(v, p) for v in x ]
end


function linear_function_fixed_offset_at_zero(x::T, p::Vector{T})::T where {T}
    return p[1] * x
end
function linear_function_fixed_offset_at_zero(x::Vector{T}, p::Vector{T})::Vector{T} where {T}
    return T[ linear_function_fixed_offset_at_zero(v, p) for v in x ]
end
