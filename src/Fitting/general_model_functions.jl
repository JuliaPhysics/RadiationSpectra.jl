FitFunction(s::Symbol) = FitFunction(Float64, Val(s))
FitFunction(T::Type{<:AbstractFloat}, s::Symbol) = FitFunction(T, Val(s))

function FitFunction(T::Type{<:AbstractFloat}, s::Val{:Gauss})
    ff = FitFunction{T}( Gauss, 1, 3 ) 
    set_parameter_names!(ff, ["A", "σ", "μ"])
    set_initial_parameters!(ff, [1, 1, 0])
    set_fitranges!(ff, ([-1, 1],))
    return ff
end
function FitFunction(T::Type{<:AbstractFloat}, s::Val{:PDF_Gauss})
    ff = FitFunction{T}( PDF_Gauss, 1, 2 ) 
    set_parameter_names!(ff, ["σ", "μ"])
    set_initial_parameters!(ff, [1, 0])
    set_fitranges!(ff, ([-1, 1],))
    return ff
end

function FitFunction(T::Type{<:AbstractFloat}, s::Val{:GaussPlusConstantBackground})
    ff = FitFunction{T}( Gauss_plus_const_background, 1, 4 ) 
    set_parameter_names!(ff, ["A", "σ", "μ", "offset"])
    set_initial_parameters!(ff, [1, 1, 0, 1])
    set_fitranges!(ff, ([-1, 1],))
    return ff
end
function FitFunction(T::Type{<:AbstractFloat}, s::Val{:GaussPlusLinearBackground})
    ff = FitFunction{T}( Gauss_plus_linear_background, 1, 5 ) 
    set_parameter_names!(ff, ["A", "σ", "μ", "offset", "lin. slope"])
    set_initial_parameters!(ff, [1, 1, 0, 1, 0])
    set_fitranges!(ff, ([-1, 1],))
    return ff
end


function Gauss(x::T, A::T, σ::T, μ::T)::T where {T <: AbstractFloat}
    r::T = A * exp( -(x - μ)^2 / (2 * σ^2)) / sqrt(2π * σ^2)
    if σ <= 0 || A < 0 return T(Inf) end
    if r < 0 r = 0 end
    return r 
end
function Gauss(x::Vector{T}, A::T, σ::T, μ::T)::Vector{T} where {T <: AbstractFloat}
    return Gauss.( x, A, σ, μ)
end
Gauss(x::Union{T, Vector{T}}, p::Vector{T}) where {T <: AbstractFloat} = Gauss(x, p... )

function PDF_Gauss(x::T, σ::T, μ::T)::T where {T <: AbstractFloat}
    r::T = exp( -(x - μ)^2 / (2 * σ^2)) / sqrt(2π * σ^2)
    if σ <= 0 return T(Inf) end
    if r < 0 r = 0 end
    return r 
end
function PDF_Gauss(x::Vector{T}, σ::T, μ::T)::Vector{T} where {T <: AbstractFloat}
    return PDF_Gauss.( x, σ, μ)
end
PDF_Gauss(x::Union{T, Vector{T}}, p::Vector{T}) where {T <: AbstractFloat} = PDF_Gauss(x, p... )

function Gauss_plus_linear_background(x::T, A::T, σ::T, μ::T, offset::T, lin_slope::T)::T where {T <: AbstractFloat}
    r::T = A * exp( -(x - μ)^2 / (2 * σ^2)) / sqrt(2π * σ^2) + lin_slope * (x - μ) + offset
    if σ <= 0 || A < 0 || offset < 0 return T(Inf) end
    if r < 0 r = 0 end
    return r 
end
function Gauss_plus_linear_background(x::Vector{T}, A::T, σ::T, μ::T, offset::T, lin_slope::T)::Vector{T} where {T <: AbstractFloat}
    return Gauss_plus_linear_background.( x, A, σ, μ, offset, lin_slope)
end
Gauss_plus_linear_background(x::Union{T, Vector{T}}, p::Vector{T}) where {T <: AbstractFloat} = Gauss_plus_linear_background(x, p... )

function Gauss_plus_const_background(x::T, A::T, σ::T, μ::T, offset::T)::T where {T <: AbstractFloat}
    r::T = A * exp( -(x - μ)^2 / (2 * σ^2)) / sqrt(2π * σ^2) + offset
    if σ <= 0 || A < 0 || offset < 0 return T(Inf) end
    if r < 0 r = 0 end
    return r
end
function Gauss_plus_const_background(x::Vector{T}, A::T, σ::T, μ::T, offset::T)::Vector{T} where {T <: AbstractFloat}
    return Gauss_plus_const_background.( x, A, σ, μ, offset)
    if σ <= 0 || A < 0 return T(Inf) end
end
Gauss_plus_const_background(x::Union{T, Vector{T}}, p::Vector{T}) where {T <: AbstractFloat} = Gauss_plus_const_background(x, p... )