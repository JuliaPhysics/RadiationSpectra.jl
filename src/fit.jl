abstract type FitBackend end
abstract type OptimFitBackend <: FitBackend end

backend_type(::Val{:optim}) = OptimFitBackend

fit(dt::Type{<:AbstractSpectrumDensity}, args...; backend::Symbol = :optim, kwargs...) = 
    fit( backend_type(Val(backend)), dt, args...; kwargs...)
fit(::Type{OptimFitBackend}, dt::Type{<:AbstractSpectrumDensity}, args...; kwargs...) = 
    opt_fit(dt, args...; kwargs...)