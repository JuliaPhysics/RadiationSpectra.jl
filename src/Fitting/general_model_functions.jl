FitFunction(s::Symbol) = FitFunction(Float64, Val(s))
FitFunction(T::Type{<:AbstractFloat}, s::Symbol) = FitFunction(T, Val(s))


function FitFunction(T::Type{<:AbstractFloat}, s::Val{:Gauss})
    ff = FitFunction{T}( Gauss, 1, 3 ) 
    set_parameter_names!(ff, ["A", "σ", "μ"])
    set_initial_parameters!(ff, [1, 1, 0])
    set_fitranges!(ff, ([-1, 1],))
    return ff
end

function FitFunction(T::Type{<:AbstractFloat}, s::Val{:GaussPlusConstantBackground})
    ff = FitFunction{T}( Gauss_plus_linear_background, 1, 4 ) 
    set_parameter_names!(ff, ["A", "σ", "μ", "offset"])
    set_initial_parameters!(ff, [1, 1, 0, 0])
    set_fitranges!(ff, ([-1, 1],))
    return ff
end
function FitFunction(T::Type{<:AbstractFloat}, s::Val{:GaussPlusLinearBackground})
    ff = FitFunction{T}( Gauss_plus_linear_background, 1, 5 ) 
    set_parameter_names!(ff, ["A", "σ", "μ", "offset", "lin. slope"])
    set_initial_parameters!(ff, [1, 1, 0, 0, 0])
    set_fitranges!(ff, ([-1, 1],))
    return ff
end


function Gauss(x::Union{T, Vector{T}}, A::T, σ::T, μ::T) where {T <: AbstractFloat}
    return @. A * exp( -(x - μ)^2 / (2 * σ^2)) / sqrt(2π * σ^2)
end
Gauss(x::Union{T, Vector{T}}, p::Vector{T}) where {T <: AbstractFloat} = Gauss(x, p... )
function PDF_Gauss(x::Union{T, Vector{T}}, σ::T, μ::T) where {T <: AbstractFloat}
    return @. exp( -(x - μ)^2 / (2 * σ^2)) / sqrt(2π * σ^2)
end
PDF_Gauss(x::Union{T, Vector{T}}, p::Vector{T}) where {T <: AbstractFloat} = PDF_Gauss(x, p... )

function Gauss_plus_linear_background(x::Union{T, Vector{T}}, A::T, σ::T, μ::T, offset::T, lin_slope::T) where {T <: AbstractFloat}
    return @. A * exp( -(x - μ)^2 / (2 * σ^2)) / sqrt(2π * σ^2) + lin_slope * (x - μ) + offset
end
Gauss_plus_linear_background(x::Union{T, Vector{T}}, p::Vector{T}) where {T <: AbstractFloat} = Gauss_plus_linear_background(x, p... )

function Gauss_plus_const_background(x::Union{T, Vector{T}}, A::T, σ::T, μ::T, offset::T) where {T <: AbstractFloat}
    return @. A * exp( -(x - μ)^2 / (2 * σ^2)) / sqrt(2π * σ^2) + offset
end
Gauss_plus_const_background(x::Union{T, Vector{T}}, p::Vector{T}) where {T <: AbstractFloat} = Gauss_plus_const_background(x, p... )