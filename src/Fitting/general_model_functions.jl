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
    set_parameter_names!(ff, ["A", "σ", "μ", "offset", "linear_slope"])
    set_initial_parameters!(ff, [1, 1, 0, 1, 0])
    set_fitranges!(ff, ([-1, 1],))
    return ff
end
function FitFunction(T::Type{<:AbstractFloat}, s::Val{:GaussPlusLinearBackgroundPlusConstBackground})
    ff = FitFunction{T}( Gauss_plus_linear_background_plus_const_background, 1, 6 ) 
    set_parameter_names!(ff, ["A", "σ", "μ", "offset", "linear_slope", "const_background"])
    set_initial_parameters!(ff, Float64[1, 1, 0, 1, 0, nextfloat(0.0)])
    set_fitranges!(ff, ([-1, 1],))
    return ff
end

function min_0(x::T) where {T <: Real}
    return x < 0 ? zero(typeof(x)) : x
end
function min_0(x)
    min_0.(x)
end

function get_inf(x::T) where {T <: Real}
    return T(Inf)
end 
function get_inf(x) 
    return get_inf.(x)
end 

function Gauss(x, p)
    scale = p[1]
    σ     = p[2]
    μ     = p[3]
    return @fastmath @. scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2))
end

function PDF_Gauss(x, p)
    σ = p[1]
    μ = p[2]
    return min_0(@fastmath @. exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2)))
end

function Gauss_plus_linear_background(x, p)
    A = p[1]
    σ = p[2]
    μ = p[3]
    offset = p[4]
    lin_slope = p[5]
    if σ <= 0 || A < 0 return get_inf(x) end
    return @fastmath @. A * exp( -(x - μ)^2 / (2 * σ^2)) / sqrt(2π * σ^2) + min_0(lin_slope * (x - μ)) + offset 
end

function Gauss_plus_linear_background_plus_const_background(x, p)
    A = p[1]
    σ = p[2]
    μ = p[3]
    offset = p[4]
    lin_slope = p[5]
    const_background = p[6]
    if σ <= 0 || A < 0 return get_inf(x) end
    return @. min_0(@fastmath @. A * exp( -(x - μ)^2 / (2 * σ^2)) / sqrt(2π * σ^2) + lin_slope * (x - μ) + offset) + const_background
end

function Gauss_plus_const_background(x, p)
    A = p[1]
    σ = p[2]
    μ = p[3]
    offset = p[4]
    if σ <= 0 || A < 0 return get_inf(x) end
    return min_0(@fastmath @. A * exp( -(x - μ)^2 / (2 * σ^2)) / sqrt(2π * σ^2) + offset)
end