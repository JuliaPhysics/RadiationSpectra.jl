abstract type AbstractSpectrumDensity{T, N} end

@inline DensityInterface.DensityKind(::AbstractSpectrumDensity) = IsDensity()

Base.eltype(d::AbstractSpectrumDensity{T}) where {T} = T

abstract type UvSpectrumDensity{T} <: AbstractSpectrumDensity{T, 1} end

function evaluate end

