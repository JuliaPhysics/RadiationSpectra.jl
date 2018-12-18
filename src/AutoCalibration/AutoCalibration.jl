function calculate_ratios(pos::Array{Float64,1})::Array{Float64, 2}
    n = length(pos)
    ratios = zeros(Float64, n, n)
    for i in 1:n
        for j in i+1:n
            ratios[i,j] = pos[i] / pos[j]
        end
    end
    return ratios
end


function determine_calibration_constant_through_peak_ratios(fPositionX::Array{<:Real,1}, photon_lines::Array{<:Real,1}; α::Float64=0.005, min_nbins::Int=50)::Tuple{Float64, Histogram}
    c::Float64 = 0
    theo_ratios = calculate_ratios(photon_lines)
    exp_ratios = calculate_ratios(fPositionX)
    n_peaks = length(fPositionX)
    n_lines = length(photon_lines)
    pcfs = Float64[] # pre calibration factors
    @inbounds for i in 1:n_peaks
        for j in i + 1:n_peaks
            for it in 1:n_lines
                for jt in it + 1:n_lines
                    if 1 - α <= theo_ratios[it,jt] / exp_ratios[i,j] <= 1 + α
                        push!(pcfs, (photon_lines[it] - photon_lines[jt]) / (fPositionX[i] - fPositionX[j]) )
                    end
                end
            end
        end
    end
    nbins = 2 * length(pcfs) > min_nbins ? 2 * length(pcfs) : min_nbins
    pcf_hist = StatsBase.fit(Histogram, pcfs, nbins= 3 * length(pcfs), closed=:left)
    c = (pcf_hist.edges[1][1:length(pcf_hist.edges[1]) - 1] .+ 0.5 * step(pcf_hist.edges[1]))[findmax(pcf_hist.weights)[2]]
    tcf = Float64[] # 'true' calibration factors
    for pcf in pcfs
        if abs(pcf-c) < step(pcf_hist.edges[1])
            push!(tcf, pcf)
        end
    end
    c = mean(tcf)
    return c, pcf_hist
end


function determine_calibration_constant_through_peak_ratios(h::Histogram{<:Real, 1, E}, photon_lines::Array{<:Real, 1}; min_n_peaks::Int=0, threshold::Real=10., α=0.005, σ::Real = 2.0)::Tuple{Float64, Histogram} where {E}
    destVector, fPositionX = peakfinder(h, threshold=threshold, sigma = σ)
    e_threshold = last(h.edges[1]) * 0.05
    delete_peak_idcs = Int[]
    for i in eachindex(fPositionX)
        if fPositionX[i] < e_threshold
            push!(delete_peak_idcs, i)
        end
    end
    deleteat!(fPositionX, delete_peak_idcs)
    while length(fPositionX) < min_n_peaks
        threshold *= 0.75
        destVector, fPositionX = peakfinder(h, threshold=threshold, sigma = σ)
        delete_peak_idcs = Int[]
        for i in eachindex(fPositionX)
            if fPositionX[i] < e_threshold
                push!(delete_peak_idcs, i)
            end
        end
        deleteat!(fPositionX, delete_peak_idcs)
    end
    c, pcf_hist= determine_calibration_constant_through_peak_ratios(fPositionX, photon_lines, α=α)
    return c, pcf_hist
end

function determine_calibration_constant_through_peak_fitting(h::Histogram{<:Real, 1, E}, photon_lines::Array{<:Real, 1}, c_pre::Real=1.0,) where {E}
    c_pre = Float64(c_pre)
    peak_fits = FitFunction[]
    for pl in photon_lines
        line::Float64 = pl 
        fitrange = (line - 20, line + 20)  
        fitrange = fitrange ./ c_pre
        first_bin = StatsBase.binindex(h, fitrange[1])
        last_bin  = StatsBase.binindex(h, fitrange[2])
        p0_sigma = 1 / c_pre 
        p0_scale = (maximum(h.weights[first_bin:last_bin]) - (h.weights[first_bin] + h.weights[last_bin]) / 2) * 2 * p0_sigma
        p0_mean = line / c_pre
        p0_bg_offset = (h.weights[first_bin] + h.weights[last_bin]) / 2
        p0_bg_slope = (h.weights[last_bin] - h.weights[first_bin]) / (fitrange[2] - fitrange[1])

        fit = FitFunction( gauss_plus_first_order_polynom )
        fit.fitrange = fitrange
        fit.initial_parameters = Float64[ p0_scale, p0_sigma, p0_mean, p0_bg_offset, p0_bg_slope ]

        lsqfit!(fit, h, estimate_uncertainties=false)
        push!(peak_fits, fit)
    end

    fitted_peak_positions = [ fr.parameters[3] for fr in peak_fits ]

    cfit = FitFunction( linear_function_fixed_offset_at_zero )
    cfit.initial_parameters = [c_pre]

    lsqfit!(cfit, fitted_peak_positions, photon_lines, estimate_uncertainties=false)
    c = cfit.parameters[1] 
    
    return c, peak_fits, cfit
end


"""
    calibrate_spectrum(h_uncal::Histogram, photon_lines::Array{Real, 1}; <keyword arguments>)::Histogram

Returns a calibrated histogram.

# Keywords
- `sigma::Real = 2.0`: The expected sigma of a peak in the spectrum. In units of bins. 
- `threshold::Real = 10.0`: Threshold for being identified as a peak in the deconvoluted spectrum. A single bin is identified as a peak when its weight exceeds the `threshold` and the previous bin was not identified as an peak.
- `min_n_peaks::Int = 0`: If the number of found peaks is smaller than `min_n_peaks` the functions lowers the parameter `threshold` until enough peaks are found.
- `α::Real = 0.005`:  = 0.5%. Acceptance level in the comparison of the peak position ratios in the peak indentification step. When the difference between the ratio of two found peak positions and the ratio of two photon lines (`photon_lines`) is smaller than `α`, the found peaks are identified as the two photon lines.

Calibrate the spectrum `h_uncal`. This is done by:
1) finding peaks through devoncolution
2) identifying them through comparison of the ratios of their positions with the ratios of the known `lines`
3) fitting all identified peaks (with a gaussian plus first order polynomial) to get their position more precisely
4) performe a linear fit (offset forced to 0) of these positions vs the true positions (`lines`) to get the calibration constant 
"""
function calibrate_spectrum(h_uncal::Histogram{<:Real, 1, E}, photon_lines::Vector{<:Real}; σ::Real = 2.0, threshold::Real = 10.0, min_n_peaks::Int = 0, α::Real = 0.005)::Histogram where {T, E}
    c_precal = determine_calibration_constant_through_peak_ratios(h_uncal, photon_lines, min_n_peaks=min_n_peaks, threshold=threshold, α=α, σ=σ)[1]
    c = determine_calibration_constant_through_peak_fitting(h_uncal, photon_lines, c_precal)[1]
    h_cal = Histogram( h_uncal.edges[1] .* c, :left )
    h_cal.weights = h_uncal.weights
    return h_cal
end 