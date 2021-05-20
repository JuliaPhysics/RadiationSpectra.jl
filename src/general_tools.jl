
_lowest_eltype(nt::NamedTuple) = _lowest_eltype(values(nt))
_lowest_eltype(nt::Tuple) = _lowest_eltype.(nt) 
_lowest_eltype(nt::T) where {T <: Real} = T
_lowest_eltype(nt::AbstractArray{T}) where {T <: Real} = T

_promote_type(nt::NamedTuple) = _promote_type(_lowest_eltype(nt)...)
_promote_type(::Type{T}, t::Tuple) where {T <: Real} = promote_type(T, t...)
_promote_type(t::Tuple, ::Type{T}) where {T <: Real} = promote_type(T, t...)
_promote_type(::Type{X}, ::Type{Y}, args... ) where {X <: Real, Y <: Real} = 
    _promote_type(promote_type(X, Y), args...)



@inline _par_in_input_form(p::RT) where {RT} = p
@inline _par_in_input_form(p::T) where {T <: ShapedAsNT} = p[]

function _remove_const_shapes_from_nt(nt::NamedTuple, shape::ValueShapes.AbstractValueShape) 
    ks = keys(nt)
    ds = filter(iks -> typeof(shape[iks]) <: ValueAccessor{<:ConstValueShape}, eachindex(ks))
    return !isempty(ds) ? delete(nt, ks[ds]) : nt 
end

function subhist(h::Histogram{<:Any, 1}, r::Tuple{<:Real,<:Real})
    first_bin, last_bin = (StatsBase.binindex(h, r[1]), StatsBase.binindex(h, r[2]))
    if first_bin < 1 first_bin = 1 end
    if (last_bin > length(h.weights)) last_bin = length(h.weights) end
    Histogram(h.edges[1][first_bin:last_bin+1], h.weights[first_bin:last_bin])
end
subhist(h, i::Interval) = subhist(h, (i.left, i.right))