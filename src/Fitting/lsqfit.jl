function lsqfit!(fit::FitFunction, xdata::Array{<:Real, 1}, ydata::Array{<:Real, 1}; estimate_uncertainties=false)::Nothing
    fr = LsqFit.curve_fit(fit.model, Float64.(xdata), Float64.(ydata), fit.initial_parameters)
    fit.parameters = fr.param
    if estimate_uncertainties fit.uncertainties = LsqFit.margin_error(fr, 1 - fit.confidence_level) end
    nothing
end
function lsqfit!(fit::FitFunction, xdata::Array{<:Real, 1}, ydata::Array{<:Real, 1}, err::Array{<:Real, 1}; estimate_uncertainties=false)::Nothing
    w = inv.(err)
    fr = LsqFit.curve_fit(fit.model, xdata, ydata, w, fit.initial_parameters)
    fit.parameters = fr.param
    if estimate_uncertainties fit.uncertainties = LsqFit.margin_error(fr, 1 - fit.confidence_level) end
    nothing
end
function lsqfit!(fit::FitFunction, xdata::Array{<:Real, 1}, ydata::Array{<:Real, 1}, xerr::Array{<:Real, 1}, yerr::Array{<:Real, 1}; estimate_uncertainties=false)::Nothing
    w = @. 1 / sqrt(xerr^2 + yerr^2)
    fr = LsqFit.curve_fit(fit.model, xdata, ydata, w, fit.initial_parameters)
    fit.parameters = fr.param
    if estimate_uncertainties fit.uncertainties = LsqFit.margin_error(fr, 1 - fit.confidence_level) end
    nothing
end

"""
    lsqfit!(fit::FitFunction, h::Histogram; estimate_uncertainties=false)::Nothing

Performs a Least Square Fit with the model `fit.model` and the initial parameters `fit.initial_parameters`
on the histogram `h` in the range `fit.fitrange`. The determined parameters are stored in `fit.parameters`
and the corresponding uncertainties in `fit.uncertainties` for the given confidence level `fit.confidence_level`.

The uncertainties are marginalizations of the covariance matrix determined by the LSQFit.jl package.
They are only calculated if the keywort `estimate_uncertainties` is set to `true`.

See [LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl) for more detail.
"""
function lsqfit!(fit::FitFunction, h::Histogram; estimate_uncertainties=false)::Nothing
    first_bin::Int = StatsBase.binindex(h, first(fit.fitrange))
    last_bin::Int  = StatsBase.binindex(h, last(fit.fitrange))
    if first_bin < 1 first_bin = 1 end
    if (last_bin > length(h.weights)) last_bin = length(h.weights) end
    x::Array{Float64,1} = collect(h.edges[1])[first_bin:last_bin] .+ 0.5 * step(h.edges[1])
    y::Array{Float64,1} = h.weights[first_bin:last_bin]
    yerr = sqrt.(y) # Poisson distributed
    yerr = [ ye != 0 ? ye : 1.  for ye in yerr] 
    fr = LsqFit.curve_fit(fit.model, x, y, yerr, fit.initial_parameters)
    fit.parameters = fr.param
    if estimate_uncertainties fit.uncertainties = LsqFit.margin_error(fr, 1 - fit.confidence_level) end
    nothing
end


