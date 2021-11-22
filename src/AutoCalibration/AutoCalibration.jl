function calculate_ratios(pos::Array{T, 1})::Array{T, 2} where {T <: AbstractFloat}
    n::Int = length(pos)
    ratios::Array{T, 2} = zeros(Float64, n, n)
    @fastmath @inbounds for i in 1:n
        for j in i+1:n
            ratios[i,j] = pos[i] / pos[j]
        end
    end
    return ratios
end
function calculate_abs_diffs(pos::Array{T, 1})::Array{T, 2} where {T <: AbstractFloat}
    n::Int = length(pos)
    ratios::Array{T, 2} = zeros(Float64, n, n)
    @fastmath @inbounds for i in 1:n
        for j in i+1:n
            ratios[i, j] = pos[j] - pos[i]
        end
    end
    return ratios
end


function determine_calibration_constant_through_peak_ratios(fPositionX::Array{<:Real,1}, photon_lines::Array{<:Real,1}; α::Float64=0.01, rtol::Real = 5e-3)
    c::Float64 = 0
    photon_lines = Float64.(photon_lines)
    photon_lines = sort(photon_lines)
    fPositionX = sort(fPositionX)
    theo_ratios = calculate_ratios(photon_lines)
    exp_ratios = calculate_ratios(fPositionX)
    Δtheo = calculate_abs_diffs(photon_lines)
    n_peaks = length(fPositionX)
    n_lines = length(photon_lines)
    pcfs = Float64[] # pre calibration factors
    calibration_constant::Float64 = 0
    Δexp::Float64 = 0
    for i in 1:n_peaks
        for j in i + 1:n_peaks
            for it in 1:n_lines-1
                for jt in it + 1:n_lines
                    if 1 - α <= theo_ratios[it, jt] / exp_ratios[i, j] <= 1 + α
                        calibration_constant = (photon_lines[it] - photon_lines[jt]) / (fPositionX[i] - fPositionX[j])
                        Δexp = calibration_constant * (fPositionX[j] - fPositionX[i])
                        exp_peak_j = fPositionX[j] * calibration_constant
                        exp_peak_i = fPositionX[i] * calibration_constant
                        peak_i_rel_diff = abs(exp_peak_i - photon_lines[it]) / photon_lines[it]
                        peak_j_rel_diff = abs(exp_peak_j - photon_lines[jt]) / photon_lines[jt]
                        if peak_i_rel_diff < rtol && peak_j_rel_diff < rtol
                            push!(pcfs, calibration_constant )
                        end
                    end
                end
            end
        end
    end
    if isempty(pcfs)
        error("No peaks identified. Try different values for keyword parameters: `α`, `rtol`")
    end

    c_abs_min_diffs = zeros(Float64, length(pcfs))

    for (ic, c) in enumerate(pcfs)
        foundlines = c * fPositionX
        for (il, l) in enumerate(photon_lines)
            diffs = (foundlines .- l) ./ l
            abs_diffs = abs.(diffs)
            min_abs_diff, nearest_idx = findmin(abs_diffs)
            c_abs_min_diffs[ic] += min_abs_diff
        end
    end
    c = pcfs[findmin(c_abs_min_diffs)[2]]
    return c
end

function _filter_peaks_from_peakfinder(h::Histogram{<:Real, 1}, peakPositions::Vector{T}, σ::Real) where {T <: Real}
    pfits = []
    for p in peakPositions
        fitrange = (p - 5σ*step(h.edges[1]), p + 5σ*step(h.edges[1]))
        d = fit(NormalPeakUvD, subhist(h, fitrange))[1]
        push!(pfits, d)
    end
    d = map(i -> abs(pfits[i].UvNormal.μ - peakPositions[i]), eachindex(peakPositions));
    σs = map(f -> f.UvNormal.σ, pfits);   
    hd = harmmean(d);
    hσ = harmmean(σs);
    accepted_peakPositions = T[]
    for (ipf, pf) in enumerate(pfits)
        if d[ipf] < 4 * hd && abs(pf.UvNormal.σ) < 2hσ
            push!(accepted_peakPositions, pf.UvNormal.μ)
        end
    end
    accepted_peakPositions
end


function determine_calibration_constant_through_peak_ratios(h::Histogram{<:Real, 1, E}, photon_lines::Array{<:Real, 1};
            min_n_peaks::Int = length(photon_lines), max_n_peaks::Int = 5 * length(photon_lines), threshold::Real = 10., α::Real = 0.01, σ::Real = 3.0, rtol::Real = 5e-3) where {E}
    h_deconv, peakPositions = peakfinder(h, threshold=threshold, σ = σ)
    fitted_peak_Positions = _filter_peaks_from_peakfinder(h, peakPositions, σ)
    while length(fitted_peak_Positions) > max_n_peaks
        threshold *= 1.1
        h_deconv, peakPositions = peakfinder(h, threshold=threshold, σ = σ)
        fitted_peak_Positions = _filter_peaks_from_peakfinder(h, peakPositions, σ)
    end
    while length(fitted_peak_Positions) < min_n_peaks
        threshold *= 0.75
        h_deconv, peakPositions = peakfinder(h, threshold=threshold, σ = σ)
        fitted_peak_Positions = _filter_peaks_from_peakfinder(h, peakPositions, σ)
    end
    c = determine_calibration_constant_through_peak_ratios(peakPositions, photon_lines, α = α, rtol = rtol)
    return c, h_deconv, peakPositions, threshold
end

function determine_calibration_constant_through_peak_fitting(h::Histogram{<:Real, 1, E}, photon_lines::Vector{<:Real}, c_pre::Real=1.0) where {E}
    c_pre = Float64(c_pre)
    peak_fits = []
    photon_lines = Float64.(photon_lines)
    processed_photon_lines = Float64[]
    for pl in photon_lines
        line::Float64 = pl
        fitrange = (line - 20, line + 20)
        fitrange = fitrange ./ c_pre

        first_bin, last_bin = StatsBase.binindex(h, fitrange[1]), StatsBase.binindex(h, fitrange[2])
        if last_bin > length(h.weights) || first_bin < 1
            @warn("$pl keV line not sufficiently contained in Histogram, skipping...");
            continue
        end
        d = fit(NormalPeakUvD, subhist(h, fitrange))[1]
        push!(peak_fits, d)
        push!(processed_photon_lines, pl)
    end
    photon_lines = intersect(photon_lines, processed_photon_lines)
    fitted_peak_positions = [ fr.UvNormal.μ for fr in peak_fits ]

    @fastmath function _linear_function_fixed_offset_at_zero(p::Vector{T}) where {T}
        return sqrt( sum(((fitted_peak_positions * p[1]) .- photon_lines).^2) )
    end
    
    opt_result = Optim.optimize( 
        _linear_function_fixed_offset_at_zero, 
        Float64[0], Float64[Inf], Float64[c_pre], 
        Fminbox(BFGS()); 
        autodiff=:forward )
    
    c = opt_result.minimizer[1]

    return c, peak_fits, opt_result
end


"""
    calibrate_spectrum(h_uncal::Histogram, photon_lines::Array{Real, 1}; <keyword arguments>)

Returns the calibrated histogram, the deconvoluted spectrum, the found (uncalibrated) peak positions and the final `threshold` value.

# Keywords
- `σ::Real = 2.0`: The expected sigma of a peak in the spectrum. In units of bins.
- `threshold::Real = 10.0`: Threshold for being identified as a peak in the deconvoluted spectrum. A single bin is identified as a peak when its weight exceeds the `threshold` and the previous bin was not identified as an peak.
- `min_n_peaks::Int = 0`: If the number of found peaks is smaller than `min_n_peaks` the functions lowers the parameter `threshold` until enough peaks are found.
- `max_n_peaks::Int = 50`: Use only the first (strongest) `max_n_peaks` peaks for peak identification.
- `α::Real = 0.005`:  = 0.5%. Acceptance level in the comparison of the peak position ratios in the peak indentification step. When the difference between the ratio of two found peak positions and the ratio of two photon lines (`photon_lines`) is smaller than `α`, the found peaks are identified as the two photon lines.
- `rtol::Real = 5e-3`:  = 5e-3. Acceptance level for tolerance of the absolute difference between true and found line position.

Calibrate the spectrum `h_uncal`. This is done by:
1) finding peaks through devoncolution
2) identifying them through comparison of the ratios of their positions with the ratios of the known `lines`
3) fitting all identified peaks (with a gaussian plus first order polynomial) to get their position more precisely
4) performe a linear fit (offset forced to 0) of these positions vs the true positions (`lines`) to get the calibration constant
"""
function calibrate_spectrum(h_uncal::Histogram{<:Real, 1}, photon_lines::Vector{<:Real}; 
            σ::Real = 3.0, threshold::Real = 50.0, min_n_peaks::Int = length(photon_lines), 
            max_n_peaks::Int = 4 * length(photon_lines), α::Real = 0.01, rtol::Real = 5e-3)
    c_precal, h_deconv, peakPositions, threshold = determine_calibration_constant_through_peak_ratios(h_uncal, photon_lines, 
                                min_n_peaks = min_n_peaks, max_n_peaks = max_n_peaks, threshold=threshold, α=α, σ=σ, rtol=rtol)
    c = determine_calibration_constant_through_peak_fitting(h_uncal, photon_lines, c_precal)[1]
    h_cal = Histogram( h_uncal.edges[1] .* c, :left )
    h_cal.weights = h_uncal.weights
    return h_cal, h_deconv, peakPositions, threshold, c, c_precal
end
