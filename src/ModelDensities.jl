abstract type AbstractModelDensity{T, N} end

eltype(d::AbstractModelDensity{T}) where {T} = T
eltype(d::Type{<:AbstractModelDensity{T}}) where {T} = T

abstract type UvModelDensity{T} <: AbstractModelDensity{T, 1} end
abstract type MvModelDensity{T,N} <: AbstractModelDensity{T, N} end

function evaluate end

