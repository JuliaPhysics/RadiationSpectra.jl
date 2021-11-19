abstract type AbstractSpectrumDistribution{T, N} end

@inline DensityInterface.DensityKind(::AbstractSpectrumDistribution) = IsDensity()

Base.eltype(d::AbstractSpectrumDistribution{T}) where {T} = T

abstract type UvSpectrumDensity{T} <: AbstractSpectrumDistribution{T, 1} end
abstract type MvSpectrumDensity{T,N} <: AbstractSpectrumDistribution{T, N} end

function evaluate end

