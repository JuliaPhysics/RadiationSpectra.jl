var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#RadiationSpectra.jl-1",
    "page": "Home",
    "title": "RadiationSpectra.jl",
    "category": "section",
    "text": "RadiationSpectra.jl provides tools to analyse Radiation Spectra.This includes (for now):Pages = [\"man/spectrum_deconvolution.md\"]Pages = [\"man/fitting.md\"]Pages = [\"man/spectrum_calibration.md\"]"
},

{
    "location": "man/spectrum_deconvolution/#",
    "page": "Peak Finding",
    "title": "Peak Finding",
    "category": "page",
    "text": ""
},

{
    "location": "man/spectrum_deconvolution/#Peak-Finding-(Spectrum-Deconvolution)-1",
    "page": "Peak Finding",
    "title": "Peak Finding (Spectrum Deconvolution)",
    "category": "section",
    "text": "See RadiationSpectra.peakfinderusing RadiationSpectra \nh_uncal = RadiationSpectra.get_example_spectrum()\nh_decon, peakpos = RadiationSpectra.peakfinder(h_uncal)\n\nusing Plots \nmyfont = Plots.font(12) # hide\npyplot(guidefont=myfont, xtickfont=myfont, ytickfont=myfont, legendfont=myfont, titlefont=myfont) # hide\np_uncal = plot(h_uncal, st=:step, label=\"Uncalibrated spectrum\", c=1, lw=1.2); \np_decon = plot(peakpos, st=:vline, c=:red, label=\"Peaks\", lw=0.3);\nplot!(h_decon, st=:step, label=\"Deconvoluted spectrum\", c=1, lw=1.2); \nplot(p_uncal, p_decon, size=(800,500), layout=(2, 1), fmt=:svg) "
},

{
    "location": "man/fitting/#",
    "page": "Fitting",
    "title": "Fitting",
    "category": "page",
    "text": ""
},

{
    "location": "man/fitting/#Fitting-1",
    "page": "Fitting",
    "title": "Fitting",
    "category": "section",
    "text": "Pages = [\"fitting.md\"]"
},

{
    "location": "man/fitting/#Likelihood-(LLH)-Fit-1D-Histogram-1",
    "page": "Fitting",
    "title": "Likelihood (LLH) Fit - 1D-Histogram",
    "category": "section",
    "text": "Get a spectrum and find a peak to fitusing Plots, RadiationSpectra, StatsBase\nmyfont = Plots.font(12) # hide\npyplot(guidefont=myfont, xtickfont=myfont, ytickfont=myfont, legendfont=myfont, titlefont=myfont) # hide\n\nh_uncal = RadiationSpectra.get_example_spectrum()\nh_decon, peakpos = RadiationSpectra.peakfinder(h_uncal)\n\nstrongest_peak_bin_idx = StatsBase.binindex(h_uncal, peakpos[1])\nstrongest_peak_bin_width = StatsBase.binvolume(h_uncal, strongest_peak_bin_idx)\nstrongest_peak_bin_amplitude = h_uncal.weights[strongest_peak_bin_idx]\n\nplot(h_uncal, st=:step, xlims=[peakpos[1] - strongest_peak_bin_width * 20, peakpos[1] + strongest_peak_bin_width * 20], size=(800,400), ylims=[0, strongest_peak_bin_amplitude * 1.1], fmt =:svg) Write a model functionfunction model(x::T, par::Vector{T}) where {T}\n    scale::T = par[1]\n    σ::T     = par[2]\n    μ::T     = par[3]\n    cp0::T   = par[4] \n    cp1::T   = par[5]\n    return scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2)) + cp0 + cp1 * (x - μ)\nend\nfunction model(x::Vector{T}, par::Vector{T}) where {T} return T[model(v, par) for v in x] end;Set up the fit function RadiationSpectra.FitFunction{T}fitfunc = RadiationSpectra.FitFunction( model ); \nfitfunc.fitrange = (peakpos[1] - 1000, peakpos[1] + 1000)\nguess_σ = strongest_peak_bin_width * 2\nguess_amplitude = strongest_peak_bin_amplitude * guess_σ\nguess_μ = peakpos[1]\nguess_offset = 0\nguess_bg_slope = 0\nfitfunc.initial_parameters = [ guess_amplitude, guess_σ, guess_μ, guess_offset, guess_bg_slope ]\nfitfuncPerforme the fit with the RadiationSpectra.llhfit!-function and plot the resultRadiationSpectra.llhfit!(fitfunc, h_uncal)\n\nplot(h_uncal, st=:step, xlims=[peakpos[1] - strongest_peak_bin_width * 20, peakpos[1] + strongest_peak_bin_width * 20], size=(800,400), label=\"Spectrum\", ylims=[0, strongest_peak_bin_amplitude * 1.1])\nplot!(fitfunc, use_initial_parameters=true, lc=:green, label=\"Guess\")\nplot!(fitfunc, lc=:red, label=\"LLH Fit\", fmt=:svg)fitfunc # hide"
},

{
    "location": "man/fitting/#LSQ-Fit-1D-Histogram-1",
    "page": "Fitting",
    "title": "LSQ Fit - 1D-Histogram",
    "category": "section",
    "text": "To perfrom a LSQ Fit on a spectrum repeat the first three steps from the Likelihood (LLH) Fit - 1D-Histogram. Then, Performe the fit with the RadiationSpectra.lsqfit!-function and plot the resultRadiationSpectra.lsqfit!(fitfunc, h_uncal)\n\nplot(h_uncal, st=:step, xlims=[peakpos[1] - strongest_peak_bin_width * 20, peakpos[1] + strongest_peak_bin_width * 20], size=(800,400), label=\"Spectrum\", ylims=[0, strongest_peak_bin_amplitude * 1.1])\nplot!(fitfunc, use_initial_parameters=true, lc=:green, label=\"Guess\")\nplot!(fitfunc, lc=:red, label=\"LSQ Fit\", fmt=:svg)fitfunc # hide"
},

{
    "location": "man/fitting/#LSQ-Fit-1D-Data-Arrays-1",
    "page": "Fitting",
    "title": "LSQ Fit - 1D-Data Arrays",
    "category": "section",
    "text": "Write a model functionusing Plots, RadiationSpectra\nmyfont = Plots.font(12) # hide\npyplot(guidefont=myfont, xtickfont=myfont, ytickfont=myfont, legendfont=myfont, titlefont=myfont) # hide\nfunction model(x::T, par::Vector{T}) where {T}\n    cp0::T   = par[1] \n    cp1::T   = par[2]\n    return cp0 + cp1 * x\nend\nfunction model(x::Vector{T}, par::Vector{T}) where {T} return T[model(v, par) for v in x] end;Create some random data xdata = Float64[1, 2, 3, 6, 8, 12]\nydata = Float64[model(x, [-0.2, 0.7]) + (rand() - 0.5) for x in xdata]\nxdata_err = Float64[0.5 for i in eachindex(xdata)]\nydata_err = Float64[1 for i in eachindex(xdata)]\n\nplot(xdata, ydata, xerr=xdata_err, yerr=ydata_err, st=:scatter, size=(800,400), label=\"Data\", fmt=:svg)Set up the fit function RadiationSpectra.FitFunction{T}fitfunc = RadiationSpectra.FitFunction( model ); \nfitfunc.fitrange = (xdata[1], xdata[end])\nguess_offset = 0\nguess_bg_slope = 1\nfitfunc.initial_parameters = Float64[ guess_offset, guess_bg_slope ]\nfitfuncPerforme the fit and plot the resultRadiationSpectra.lsqfit!(fitfunc, xdata, ydata, xdata_err, ydata_err) # xdata_err and ydata_err are optional\n\nplot(xdata, ydata, xerr=xdata_err, yerr=ydata_err, st=:scatter, size=(800,400), label=\"Data\")\nplot!(fitfunc, use_initial_parameters=true, lc=:green, label=\"Guess\")\nplot!(fitfunc, lc=:red, label=\"LSQ Fit\", fmt=:svg)fitfunc # hide"
},

{
    "location": "man/spectrum_calibration/#",
    "page": "Calibration",
    "title": "Calibration",
    "category": "page",
    "text": ""
},

{
    "location": "man/spectrum_calibration/#Spectrum-Calibration-1",
    "page": "Calibration",
    "title": "Spectrum Calibration",
    "category": "section",
    "text": "The function RadiationSpectra.calibrate_spectrum return a calibrated spectrum of an uncalibrated one.  As input it needs the uncalibrated spectrum (StatsBase::Histogram) and a Vector of known peak positions like photon lines.The calibration is based on the assumption that the calibration function is just f(x) = ccdot x.using Plots, RadiationSpectra \nmyfont = Plots.font(12) # hide\npyplot(guidefont=myfont, xtickfont=myfont, ytickfont=myfont, legendfont=myfont, titlefont=myfont) # hide\n\nh_uncal = RadiationSpectra.get_example_spectrum()\nphoton_lines = [609.312, 911.204, 1120.287, 1460.830, 1764.494] # keV\nh_cal = RadiationSpectra.calibrate_spectrum(h_uncal, photon_lines)\n\np_uncal = plot(h_uncal, st=:step, label=\"Uncalibrated spectrum\"); \np_cal = plot(h_cal, st=:step, label=\"Calibrated spectrum\", xlabel=\"E / keV\"); \nvline!(photon_lines, lw=0.5, color=:red, label=\"Photon lines\")\nplot(p_uncal, p_cal, size=(800,500), layout=(2, 1), fmt=:svg) Zoom into one of the peaks:p_cal = plot(h_cal, st=:step, label=\"Calibrated spectrum\", xlabel=\"E / keV\", xlims=[1440, 1480], size=[800, 400]); # hide\nvline!([1460.830], label=\"Photon line\", fmt=:svg) # hide"
},

{
    "location": "man/spectrum_calibration/#Algorithm-1",
    "page": "Calibration",
    "title": "Algorithm",
    "category": "section",
    "text": "Deconvolution -> Peak finding \nPeak Identification - Which peak corresponds to which photon line? \nFit identified peaks \nFit determined peak positions vs \'true\' positions (photon lines) to get the calibration constant c."
},

{
    "location": "api/#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "api/#API-1",
    "page": "API",
    "title": "API",
    "category": "section",
    "text": "DocTestSetup  = quote\n    using RadiationSpectra\nend"
},

{
    "location": "api/#Types-1",
    "page": "API",
    "title": "Types",
    "category": "section",
    "text": "Order = [:type]"
},

{
    "location": "api/#Functions-1",
    "page": "API",
    "title": "Functions",
    "category": "section",
    "text": "Order = [:function]"
},

{
    "location": "api/#RadiationSpectra.peakfinder-Tuple{StatsBase.Histogram}",
    "page": "API",
    "title": "RadiationSpectra.peakfinder",
    "category": "method",
    "text": "peakfinder(h::Histogram; <keyword arguments>)::Tuple{Histogram, Array{Float64, 1}}\n\nReturns a deconvoluted spectrum and an array of peak positions.\n\nKeywords\n\nsigma::Real=2.0: The expected sigma of a peak in the spectrum. In units of bins. \nthreshold::Real=10.0: Threshold for being identified as a peak in the deconvoluted spectrum. A single bin is identified as a peak when its weight exceeds the threshold and the previous bin was not identified as an peak.\nbackgroundRemove::Bool=true\ndeconIterations::Int=3\nmarkov::Bool=true\naverWindow::Int=3\n\nSource\n\nThis function is basically a copy of TSpectrum::SearchHighRes from ROOT.\n\nM.A. Mariscotti: A method for identification of peaks in the presence of background and its application to spectrum analysis. NIM 50 (1967), 309-320.\nM. Morhac;, J. Kliman, V. Matouoek, M. Veselsky, I. Turzo.:Identification of peaks in multidimensional coincidence gamma-ray spectra. NIM, A443 (2000) 108-125.\nZ.K. Silagadze, A new algorithm for automatic photopeak searches. NIM A 376 (1996), 451.\n\n\n\n\n\n"
},

{
    "location": "api/#RadiationSpectra.AbstractFitFunction",
    "page": "API",
    "title": "RadiationSpectra.AbstractFitFunction",
    "category": "type",
    "text": "AbstractFitFunction{T, N}\n\nAbstract type for an N-dimensional fit of eltype T.\n\n\n\n\n\n"
},

{
    "location": "api/#RadiationSpectra.FitFunction",
    "page": "API",
    "title": "RadiationSpectra.FitFunction",
    "category": "type",
    "text": "FitFunction{T} <: AbstractFitFunction{T, 1}\n\nFields:\n\nmodel::Function: Function of the fit model.\nfitrange::Tuple{Union{Missing, T}, Union{Missing, T}}: Range on which the fit is performed.\nparameters::Union{Vector{Missing}, Vector{T}}: Fitted parameters.\nuncertainties::Union{Vector{Missing}, Vector{T}}: Uncertainties of the fitted parameters.\ninitial_parameters::Union{Vector{Missing}, Vector{T}}: Initial parameters.\nconfidence_level::Union{Missing, T}: Confidence level, used to determine the uncertainties.\n\nPlotting:\n\nA plot recipes exists for this struct: plot(fit::FitFunction{T}; npoints=501)\n\nPlots the model function over the fitrange with 501-points. \n\n\n\n\n\n"
},

{
    "location": "api/#RadiationSpectra.calibrate_spectrum-Union{Tuple{E}, Tuple{T}, Tuple{Histogram{#s118,1,E} where #s118<:Real,Array{#s117,1} where #s117<:Real}} where E where T",
    "page": "API",
    "title": "RadiationSpectra.calibrate_spectrum",
    "category": "method",
    "text": "calibrate_spectrum(h_uncal::Histogram, photon_lines::Array{Real, 1}; <keyword arguments>)::Histogram\n\nReturns a calibrated histogram.\n\nKeywords\n\nsigma::Real = 2.0: The expected sigma of a peak in the spectrum. In units of bins. \nthreshold::Real = 10.0: Threshold for being identified as a peak in the deconvoluted spectrum. A single bin is identified as a peak when its weight exceeds the threshold and the previous bin was not identified as an peak.\nmin_n_peaks::Int = 0: If the number of found peaks is smaller than min_n_peaks the functions lowers the parameter threshold until enough peaks are found.\nα::Real = 0.005:  = 0.5%. Acceptance level in the comparison of the peak position ratios in the peak indentification step. When the difference between the ratio of two found peak positions and the ratio of two photon lines (photon_lines) is smaller than α, the found peaks are identified as the two photon lines.\n\nCalibrate the spectrum h_uncal. This is done by:\n\nfinding peaks through devoncolution\nidentifying them through comparison of the ratios of their positions with the ratios of the known lines\nfitting all identified peaks (with a gaussian plus first order polynomial) to get their position more precisely\nperforme a linear fit (offset forced to 0) of these positions vs the true positions (lines) to get the calibration constant \n\n\n\n\n\n"
},

{
    "location": "api/#RadiationSpectra.gauss-Union{Tuple{T}, Tuple{Array{T,1},Array{T,1}}} where T",
    "page": "API",
    "title": "RadiationSpectra.gauss",
    "category": "method",
    "text": "gauss(x::Vector{T}, p::Vector{T})::Vector{T}\n\nMaps x to gauss(x::T, p::Vector{T})::T\n\n\n\n\n\n"
},

{
    "location": "api/#RadiationSpectra.gauss-Union{Tuple{T}, Tuple{T,Array{T,1}}} where T",
    "page": "API",
    "title": "RadiationSpectra.gauss",
    "category": "method",
    "text": "gauss(x::T, p::Vector{T})::T\n\nA Gaussian with 3 parameters:\n\np[1]: Scale/Amplitude\np[2]: σ\np[3]: μ\n\n\n\n\n\n"
},

{
    "location": "api/#RadiationSpectra.get_example_spectrum-Tuple{}",
    "page": "API",
    "title": "RadiationSpectra.get_example_spectrum",
    "category": "method",
    "text": "get_example_spectrum()::Histogram\n\nReturns an uncalibrated radiation spectrum for testing and demonstrating purpose.\n\n\n\n\n\n"
},

{
    "location": "api/#RadiationSpectra.llhfit!-Tuple{RadiationSpectra.FitFunction,StatsBase.Histogram}",
    "page": "API",
    "title": "RadiationSpectra.llhfit!",
    "category": "method",
    "text": " llhfit!(fit::FitFunction, h::Histogram)::Nothing\n\nPerforms a log-likelihood fit of the model function fit.model and the initial parameters fit.initial_parameters on the histogram h in the range fit.fitrange. The determined parameters are stored in fit.parameters.\n\nThe likelihood for each individual bin is the Poission distribution.\n\nThere are no uncertainty estimations for this fit yet.\n\n\n\n\n\n"
},

{
    "location": "api/#RadiationSpectra.lsqfit!-Tuple{RadiationSpectra.FitFunction,StatsBase.Histogram}",
    "page": "API",
    "title": "RadiationSpectra.lsqfit!",
    "category": "method",
    "text": "lsqfit!(fit::FitFunction, h::Histogram; estimate_uncertainties=false)::Nothing\n\nPerforms a Least Square Fit with the model fit.model and the initial parameters fit.initial_parameters on the histogram h in the range fit.fitrange. The determined parameters are stored in fit.parameters and the corresponding uncertainties in fit.uncertainties for the given confidence level fit.confidence_level.\n\nThe uncertainties are marginalizations of the covariance matrix determined by the LSQFit.jl package. They are only calculated if the keywort estimate_uncertainties is set to true.\n\nSee LsqFit.jl for more detail.\n\n\n\n\n\n"
},

{
    "location": "api/#Documentation-1",
    "page": "API",
    "title": "Documentation",
    "category": "section",
    "text": "Modules = [RadiationSpectra]\nOrder = [:type, :function]"
},

{
    "location": "LICENSE/#",
    "page": "LICENSE",
    "title": "LICENSE",
    "category": "page",
    "text": ""
},

{
    "location": "LICENSE/#LICENSE-1",
    "page": "LICENSE",
    "title": "LICENSE",
    "category": "section",
    "text": "using Markdown\nMarkdown.parse_file(joinpath(@__DIR__, \"..\", \"..\", \"LICENSE.md\"))"
},

]}
