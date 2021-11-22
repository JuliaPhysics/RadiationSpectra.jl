abstract type FitBackend end
abstract type OptimFitBackend <: FitBackend end

backend_type(::Val{:optim}) = OptimFitBackend

fit(args...; backend::Symbol = :optim, kwargs...) = fit( backend_type(Val(backend)), args...; kwargs...)
fit(::Type{OptimFitBackend}, args...; kwargs...) = opt_fit(args...; kwargs...)