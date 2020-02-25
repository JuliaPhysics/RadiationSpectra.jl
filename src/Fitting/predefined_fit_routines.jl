function fit_single_peak_histogram(h::StatsBase.Histogram, hist_data::Union{Vector{T}, Missing} = missing; fit_function = :Gauss, weights = missing ) where T<:Real
    bins = collect(h.edges[1])[1:end-1]
    max = maximum( h.weights)
    binsize = bins[2]-bins[1]
    bins .+= 0.5 * binsize
    mymean = bins[findfirst(x->x==maximum(h.weights), h.weights)]
    sigma = ismissing(hist_data) ? 3 * binsize : stdm(filter(x-> isapprox(x, mymean, atol = sqrt(abs(mymean))), hist_data), mymean)
    ampl = max / abs(sigma)
    p0 = [ampl, sigma, mymean]
    fit_obj = FitFunction(fit_function)
    set_fitranges!(fit_obj,([bins[1], bins[end]],))
    set_initial_parameters!(fit_obj, p0)
    lsqfit!(fit_obj, h, weights = weights)
    fit_obj
end

function fit_single_peak_histogram_w_gaussian(h::StatsBase.Histogram, hist_data::Union{Vector{T}, Missing} = missing; weights = missing) where T<:Real
    fit_single_peak_histogram(h, hist_data, fit_function = :Gauss, weights = weights)
end

function fit_single_peak_histogram_w_cauchy(h::StatsBase.Histogram, hist_data::Union{Vector{T}, Missing} = missing; weights = missing) where T<:Real
    fit_single_peak_histogram(h, hist_data, fit_function = :Cauchy, weights = weights)
end

function fit_single_peak_histogram_refined(h::StatsBase.Histogram, hist_data::Union{Vector{T}, Missing} = missing; fit_function::Symbol = :Gauss, weights = missing ) where T<:Real
    @assert fit_function in [:Gauss,:Cauchy] "Please specify either :Gauss or :Cauchy. Initial Parameter Estimation is attempted for you."
    bins = collect(h.edges[1])[1:end-1]
    max = maximum( h.weights)
    binsize = bins[2]-bins[1]
    bins .+= 0.5 * binsize
    mymean = bins[findfirst(x->x==maximum(h.weights), h.weights)]
    sigma = ismissing(hist_data) ? 3 * binsize : stdm(filter(x-> isapprox(x, mymean, atol = sqrt(abs(mymean))), hist_data), mymean)
    ampl = max / abs(sigma)
    p0 = [ampl, sigma, mymean]
    fit_obj = FitFunction(fit_function)
    set_fitranges!(fit_obj,([bins[1], bins[end]],))
    set_initial_parameters!(fit_obj, p0)
    lsqfit!(fit_obj, h)
    set_fitranges!(fit_obj,([fit_obj.fitted_parameters[3] - 6*abs(fit_obj.fitted_parameters[2]), fit_obj.fitted_parameters[3] + 6*abs(fit_obj.fitted_parameters[2])],))
    set_initial_parameters!(fit_obj, fit_obj.fitted_parameters)
    lsqfit!(fit_obj, h, weights = weights)
    fit_obj
end

function fit_single_peak_histogram_w_gaussian_refined(h::StatsBase.Histogram, hist_data::Union{Vector{T}, Missing} = missing; weights = missing) where T<:Real
    fit_single_peak_histogram_refined(h, hist_data, fit_function = :Gauss, weights = weights)
end

function fit_single_peak_histogram_w_cauchy_refined(h::StatsBase.Histogram, hist_data::Union{Vector{T}, Missing} = missing; weights = missing) where T<:Real
    fit_single_peak_histogram_refined(h, hist_data, fit_function = :Cauchy, weights = weights)
end
