"""
    fit_single_peak_histogram(h::StatsBase.Histogram, fitrange::Union{Tuple{T,T}, Missing} = missing; fit_function = :Gauss, weights = missing )::FitFunction{T, 1, NP} where {T <: AbstractFloat, NP}

This method is meant for quick and stable fitting. The estimation of initial parameters is done automatically, this allows and easy fit with just one line of code. The user is expected to pass a histogram containing only one peak or specify corresponding fitrange.
The method performs a least square fit with the model specified under  `fit_function = ...`. Currently supported are :Gauss (default), :Gauss_pol1, Cauchy.
"""
function fit_single_peak_histogram(h::StatsBase.Histogram, fitrange::Union{Tuple{T,T}, Missing} = missing; fit_function = :Gauss, weights = missing ) where T<:Real
    first_bin, last_bin = ismissing(fitrange) ? (1, length(h.weights)) : (StatsBase.binindex(h, fitrange[1]), StatsBase.binindex(h, fitrange[2]))
    if first_bin < 1 first_bin = 1 end
    if (last_bin > length(h.weights)) last_bin = length(h.weights) end
    hist_weights = h.weights[first_bin:last_bin]
    bin_centers = StatsBase.midpoints(h.edges[1])[first_bin:last_bin]
    bin_widths = map(i -> StatsBase.binvolume(h, i), first_bin:last_bin)

    p0_mean, p0_sigma = mean_and_std(bin_centers, FrequencyWeights(hist_weights), corrected=true)
    p0_scale = sum(h.weights[StatsBase.binindex(h, p0_mean-p0_sigma):StatsBase.binindex(h, p0_mean+p0_sigma)]) / 0.68

    p0, fit_obj = if fit_function in [:Gauss, :Cauchy]
            [p0_scale, p0_sigma, p0_mean], FitFunction(fit_function)
        elseif fit_function in [:Gauss_pol1]
            p0_bg_offset = (hist_weights[1]/bin_widths[1] + hist_weights[end]/bin_widths[end]) / 2
            p0_bg_slope = (hist_weights[1]/bin_widths[1] - hist_weights[end]/bin_widths[end]) / (fitrange[2] - fitrange[1])
            [p0_scale, p0_sigma, p0_mean, p0_bg_offset, p0_bg_slope], FitFunction(:GaussPlusLinearBackground)
    end
    set_fitranges!(fit_obj,([bin_centers[1], bin_centers[end]],))
    set_initial_parameters!(fit_obj, p0)
    lsqfit!(fit_obj, h, weights = weights)
    fit_obj
end


function fit_single_peak_histogram_w_gaussian(h::StatsBase.Histogram, fitrange = missing; weights = missing) where T<:Real
    fit_single_peak_histogram(h, fitrange, fit_function = :Gauss, weights = weights)
end

function fit_single_peak_histogram_w_cauchy(h::StatsBase.Histogram, fitrange = missing; weights = missing) where T<:Real
    fit_single_peak_histogram(h, fitrange, fit_function = :Cauchy, weights = weights)
end

"""
    fit_single_peak_histogram_refined(h::StatsBase.Histogram, fitrange::Union{Tuple{T,T}, Missing} = missing; n_sig = 3.0, fit_function = :Gauss, weights = missing )::FitFunction{T, 1, NP} where {T <: AbstractFloat, NP}

This method is meant for quick and stable fitting. The estimation of initial parameters is done automatically, this allows and easy fit with just one line of code. The user is expected to pass a histogram containing only one peak or specify corresponding fitrange.
The method performs a least square fit with the model specified under  `fit_function = ...`. Currently supported are :Gauss (default), :Gauss_pol1, Cauchy. It then fits again using the results of the first fit as initial parameters.
The fitrange for the second fit can be customized using the `n_sig = ...` keyword.
"""
function fit_single_peak_histogram_refined(h::StatsBase.Histogram, fitrange::Union{Tuple{T,T}, Missing} = missing; n_sig = 4.0, fit_function::Symbol = :Gauss, weights = missing ) where T<:Real
    @assert fit_function in [:Gauss, :Gauss_pol1, :Cauchy] "Please specify 'fit_function = ' either :Gauss, :Gauss_pol1 or :Cauchy as keyword argument. Initial Parameter Estimation is attempted for you."
    fit_obj = fit_single_peak_histogram(h, fitrange, fit_function = fit_function, weights = weights)
    n_sig = minimum([abs(fit_obj.fitted_parameters[3] - fit_obj.bin_centers[1]) / abs(fit_obj.fitted_parameters[2]) , abs(fit_obj.bin_centers[end] - fit_obj.fitted_parameters[3]) / abs(fit_obj.fitted_parameters[2]), n_sig])
    set_fitranges!(fit_obj,([fit_obj.fitted_parameters[3] - n_sig*abs(fit_obj.fitted_parameters[2]), fit_obj.fitted_parameters[3] + n_sig*abs(fit_obj.fitted_parameters[2])],))
    set_initial_parameters!(fit_obj, fit_obj.fitted_parameters)
    lsqfit!(fit_obj, h, weights = weights)
    fit_obj
end

function fit_single_peak_histogram_w_gaussian_refined(h::StatsBase.Histogram, fitrange = missing; weights = missing) where T<:Real
    fit_single_peak_histogram_refined(h, fitrange, fit_function = :Gauss, weights = weights)
end

function fit_single_peak_histogram_w_cauchy_refined(h::StatsBase.Histogram, fitrange = missing; weights = missing) where T<:Real
    fit_single_peak_histogram_refined(h, fitrange, fit_function = :Cauchy, weights = weights)
end
