var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"DocTestSetup  = quote\n    using RadiationSpectra\nend","category":"page"},{"location":"api/#Types","page":"API","title":"Types","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Order = [:type]","category":"page"},{"location":"api/#Functions","page":"API","title":"Functions","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Order = [:function]","category":"page"},{"location":"api/#Documentation","page":"API","title":"Documentation","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [RadiationSpectra]\nOrder = [:type, :function]","category":"page"},{"location":"api/#RadiationSpectra.calibrate_spectrum-Tuple{StatsBase.Histogram{<:Real, 1}, Vector{<:Real}}","page":"API","title":"RadiationSpectra.calibrate_spectrum","text":"calibrate_spectrum(h_uncal::Histogram, photon_lines::Array{Real, 1}; <keyword arguments>)\n\nReturns the calibrated histogram, the deconvoluted spectrum, the found (uncalibrated) peak positions and the final threshold value.\n\nKeywords\n\nσ::Real = 2.0: The expected sigma of a peak in the spectrum. In units of bins.\nthreshold::Real = 10.0: Threshold for being identified as a peak in the deconvoluted spectrum. A single bin is identified as a peak when its weight exceeds the threshold and the previous bin was not identified as an peak.\nmin_n_peaks::Int = 0: If the number of found peaks is smaller than min_n_peaks the functions lowers the parameter threshold until enough peaks are found.\nmax_n_peaks::Int = 50: Use only the first (strongest) max_n_peaks peaks for peak identification.\nα::Real = 0.005:  = 0.5%. Acceptance level in the comparison of the peak position ratios in the peak indentification step. When the difference between the ratio of two found peak positions and the ratio of two photon lines (photon_lines) is smaller than α, the found peaks are identified as the two photon lines.\nrtol::Real = 5e-3:  = 5e-3. Acceptance level for tolerance of the absolute difference between true and found line position.\n\nCalibrate the spectrum h_uncal. This is done by:\n\nfinding peaks through devoncolution\nidentifying them through comparison of the ratios of their positions with the ratios of the known lines\nfitting all identified peaks (with a gaussian plus first order polynomial) to get their position more precisely\nperforme a linear fit (offset forced to 0) of these positions vs the true positions (lines) to get the calibration constant\n\n\n\n\n\n","category":"method"},{"location":"api/#RadiationSpectra.get_example_spectrum-Tuple{}","page":"API","title":"RadiationSpectra.get_example_spectrum","text":"get_example_spectrum()::Histogram\n\nReturns an uncalibrated radiation spectrum for testing and demonstrating purpose.\n\n\n\n\n\n","category":"method"},{"location":"api/#RadiationSpectra.opt_fit","page":"API","title":"RadiationSpectra.opt_fit","text":"# opt_fit(DT::Type{<:UvSpectrumDensity}, h::Histogram{<:Any, 1},\n         p0::AbstractVector, lower_bounds::AbstractVector, upper_bounds::AbstractVector)\n\nMaximum Likelihood Estimation Fit of model density d on the histogram h.     \n\n\n\n\n\n","category":"function"},{"location":"api/#RadiationSpectra.opt_fit-2","page":"API","title":"RadiationSpectra.opt_fit","text":"# opt_fit(DT::Type{<:UvSpectrumDensity}, h::Histogram{<:Any, 1},\n         p0::NamedTuple, lower_bounds::NamedTuple, upper_bounds::NamedTuple)\n\nMaximum Likelihood Estimation Fit of model density d on the histogram h.     \n\n\n\n\n\n","category":"function"},{"location":"api/#RadiationSpectra.peakfinder-Tuple{StatsBase.Histogram}","page":"API","title":"RadiationSpectra.peakfinder","text":"peakfinder(h::Histogram; <keyword arguments>)::Tuple{Histogram, Array{Float64, 1}}\n\nReturns a deconvoluted spectrum and an array of peak positions.\n\nKeywords\n\nσ::Real=2.0: The expected sigma of a peak in the spectrum. In units of bins. \nthreshold::Real=10.0: Threshold for being identified as a peak in the deconvoluted spectrum. A single bin is identified as a peak when its weight exceeds the threshold and the previous bin was not identified as an peak.\nbackgroundRemove::Bool=true\ndeconIterations::Int=3\nmarkov::Bool=true\naverWindow::Int=3\n\nSource\n\nThis function is basically a copy of TSpectrum::SearchHighRes from ROOT.\n\nM.A. Mariscotti: A method for identification of peaks in the presence of background and its application to spectrum analysis. NIM 50 (1967), 309-320.\nM. Morhac;, J. Kliman, V. Matouoek, M. Veselsky, I. Turzo.:Identification of peaks in multidimensional coincidence gamma-ray spectra. NIM, A443 (2000) 108-125.\nZ.K. Silagadze, A new algorithm for automatic photopeak searches. NIM A 376 (1996), 451.\n\n\n\n\n\n","category":"method"},{"location":"man/spectrum_calibration/#Spectrum-Calibration","page":"Calibration","title":"Spectrum Calibration","text":"","category":"section"},{"location":"man/spectrum_calibration/","page":"Calibration","title":"Calibration","text":"The function RadiationSpectra.calibrate_spectrum return a calibrated spectrum of an uncalibrated one.  As input it needs the uncalibrated spectrum (StatsBase::Histogram) and a Vector of known peak positions like photon lines.","category":"page"},{"location":"man/spectrum_calibration/","page":"Calibration","title":"Calibration","text":"The calibration is based on the assumption that the calibration function is just f(x) = ccdot x.","category":"page"},{"location":"man/spectrum_calibration/","page":"Calibration","title":"Calibration","text":"using Plots, RadiationSpectra \ngr(lw = 2.0); # hide\n\nh_uncal = RadiationSpectra.get_example_spectrum()\nphoton_lines = [609.312, 911.204, 1120.287, 1460.830, 1764.494] # keV\nh_cal, h_deconv, peakPositions, threshold, c, c_precal = RadiationSpectra.calibrate_spectrum(h_uncal, photon_lines)\n\np_uncal = plot(h_uncal, st=:step, label=\"Uncalibrated spectrum\"); \np_deconv = plot(h_deconv, st=:step, label = \"Deconvoluted spectrum\");\nhline!([threshold], label = \"threshold\", lw = 1.5);\np_cal = plot(h_cal, st=:step, label=\"Calibrated spectrum\", xlabel=\"E / keV\"); \nvline!(photon_lines, lw=0.5, color=:red, label=\"Photon lines\");\nplot(p_uncal, p_deconv, p_cal, size=(800,700), layout=(3, 1));\nsavefig(\"calibration_of_spectrum.svg\"); # hide\nsavefig(\"calibration_of_spectrum.pdf\"); nothing # hide","category":"page"},{"location":"man/spectrum_calibration/","page":"Calibration","title":"Calibration","text":"(Image: Data)","category":"page"},{"location":"man/spectrum_calibration/","page":"Calibration","title":"Calibration","text":"Zoom into one of the peaks:","category":"page"},{"location":"man/spectrum_calibration/","page":"Calibration","title":"Calibration","text":"p_cal = plot(h_cal, st=:step, label=\"Calibrated spectrum\", xlabel=\"E / keV\", xlims=[1440, 1480], size=[800, 400]); # hide\nvline!([1460.830], label=\"Photon line\"); # hide\nsavefig(\"calibration_of_spectrum_zoom.svg\"); # hide\nsavefig(\"calibration_of_spectrum_zoom.pdf\"); nothing # hide","category":"page"},{"location":"man/spectrum_calibration/","page":"Calibration","title":"Calibration","text":"(Image: Data)","category":"page"},{"location":"man/spectrum_calibration/#Algorithm","page":"Calibration","title":"Algorithm","text":"","category":"section"},{"location":"man/spectrum_calibration/","page":"Calibration","title":"Calibration","text":"Deconvolution -> Peak finding \nPeak Identification - Which peak corresponds to which photon line? \nFit identified peaks \nFit determined peak positions vs 'true' positions (photon lines) to get the calibration constant c.","category":"page"},{"location":"man/fitting/#Fitting","page":"Fitting","title":"Fitting","text":"","category":"section"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"Pages = [\"fitting.md\"]","category":"page"},{"location":"man/fitting/#Likelihood-(LLH)-Fit-1D-Histogram","page":"Fitting","title":"Likelihood (LLH) Fit - 1D-Histogram","text":"","category":"section"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"The underlying distribution for the Likelihood of a model to a histogram  is the Poisson distribution. The evaluation of the likelihood of an model with a set of parameters on the data in form of a 1D-Histogram is implemented in this package. ","category":"page"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"Here, it is shown exemplary, how a model is defined and how it is fitted to a histogram. ","category":"page"},{"location":"man/fitting/#.-The-Data","page":"Fitting","title":"1. The Data","text":"","category":"section"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"The data will be in form of 1D-Histograms.  In this example, the histogram is a small subhistogram around a peak in the uncalibrated spectrum of a germanium detector. ","category":"page"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"using Plots, RadiationSpectra, StatsBase, Distributions\ngr(lw = 2.0); # hide\n\nh_uncal = RadiationSpectra.get_example_spectrum()\nh_decon, peakpos = RadiationSpectra.peakfinder(h_uncal)\n\nstrongest_peak_bin_idx = StatsBase.binindex(h_uncal, peakpos[1])\nstrongest_peak_bin_width = StatsBase.binvolume(h_uncal, strongest_peak_bin_idx)\nstrongest_peak_bin_amplitude = h_uncal.weights[strongest_peak_bin_idx]\n\nh_sub = RadiationSpectra.subhist(h_uncal, (peakpos[1] - strongest_peak_bin_width * 20, peakpos[1] + strongest_peak_bin_width * 20))\n\nplot(h_sub, st=:step, size=(800,400));\nsavefig(\"only_data_hist.pdf\"); # hide\nsavefig(\"only_data_hist.svg\"); nothing # hide","category":"page"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"(Image: Data)","category":"page"},{"location":"man/fitting/#.-User-Defined-Spectrum-Density","page":"Fitting","title":"2. User Defined Spectrum Density","text":"","category":"section"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"Now we want to define the density, which we want to fit to the data. ","category":"page"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"It is a new struct which needs to be subtype of RadiationSpectra.UvSpectrumDensity{T}. Here, our model will be a Gaussian (signal) on top of flat offset (background):","category":"page"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"struct CustomSpectrumDensity{T} <: RadiationSpectra.UvSpectrumDensity{T}\n    A::T\n    σ::T\n    μ::T\n    offset::T\nend","category":"page"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"Also, a method for the this type for the function RadiationSpectra.evaluate(d::UvSpectrumDensity, x) has to be defined. Note, that these are supposed to be densities as the returned value  are multiplied to the bin volume of the corresponding bin in the calculation of the likelihood. ","category":"page"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"function RadiationSpectra.evaluate(d::CustomSpectrumDensity, x)\n    return d.A * pdf(Normal(d.μ, d.σ), x) + d.offset\nend","category":"page"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"The final step is to define a constructor for the model for a set of parameters.  Here, the parameters are in form of NamedTuple, which is the recommended form. But, also a simple vector can be used.  The package ValueShapes.jl is integrated into RadiationSpectra.jl.  Also vectors or matrices can be used within NamedTuple's. Here, only scalars will be used: ","category":"page"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"function CustomSpectrumDensity(nt::NamedTuple{(:A,:μ,:σ,:offset)})\n    T = promote_type(typeof.(values(nt))...)\n    CustomSpectrumDensity(T(nt.A), T(nt.σ), T(nt.μ), T(nt.offset))\nend\nnothing # hide","category":"page"},{"location":"man/fitting/#.-Set-up-Initial-Guess-and-Bounds-of-the-Parameters","page":"Fitting","title":"3. Set up Initial Guess and Bounds of the Parameters","text":"","category":"section"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"Either a maximum likelihood estimation (MLE) (default) or and bayesian fit can be performed.  In the maximum likelihood fit the package Optim.jl is used to maximize the likelihood.  The bayesian fit is performed via the BAT.jl package.","category":"page"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"In order to perform the MLE an initial guess, lower bounds and upper bounds for the parameter have to be given:","category":"page"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"p0 = (A = 4000.0, μ = 129500, σ = 500, offset = 100) \nlower_bounds = (A = 0.0, μ = 128000, σ = 0.1, offset = 0) \nupper_bounds = (A = 30000.0, μ = 131000, σ = 1000, offset = 500) \nnothing # hide","category":"page"},{"location":"man/fitting/#.-Perform-the-Fit","page":"Fitting","title":"4. Perform the Fit","text":"","category":"section"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"Now we are ready to perform the fit.  The syntax follows the syntax of the fit function of  StatsBase.jl and  Distributions.jl: fit(::Type{Model}, data)::Model. Here, the Model will by our just defined model CustomSpectrumDensity and the data will be the histogram h_sub. Also the initial guess and bounds have to be parsed as additional arguments:","category":"page"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"fitted_dens, backend_result = fit(CustomSpectrumDensity, h_sub, p0, lower_bounds, upper_bounds)","category":"page"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"The first returned argument, fitted_dens, is an instance of CustomSpectrumDensity with the fitted parameters. ","category":"page"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"The second returned argument, backend_result, is, in case of a MLE fit, the returned object of the optimizer  of Optim.jl used to perform the maximization of the likelihood.  In case of a bayesian fit, via BAT.jl as fit backend, the second argument would be samples from BAT.jl.","category":"page"},{"location":"man/fitting/#.-Plot-the-fitted-Density:","page":"Fitting","title":"5. Plot the fitted Density:","text":"","category":"section"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"The package provides a simple plot recipe to plot the fitted density on top of the data through the defined evaluate method and the binning of the histogram:","category":"page"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"plot(h_sub, st=:step, size=(800,400), label=\"Spectrum\");\nplot!(fitted_dens, h_sub, label = \"Fit\");\nsavefig(\"data_hist_plus_fit.pdf\"); # hide\nsavefig(\"data_hist_plus_fit.svg\"); nothing # hide","category":"page"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"(Image: Data)","category":"page"},{"location":"man/fitting/#Optional:-Fixen-Parameters-in-the-Fit","page":"Fitting","title":"Optional: Fixen Parameters in the Fit","text":"","category":"section"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"The package ValueShapes.jl can be used to set individual parameters constant in the fit by defining the shape of parameters. E.g.:","category":"page"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"using ValueShapes\nshape = NamedTupleShape(\n    A = ScalarShape{Real}(),\n    μ = ScalarShape{Real}(),\n    σ = ScalarShape{Real}(),\n    offset = ConstValueShape{Real}(p0.offset),\n)","category":"page"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"The shape has to be passed to the fit function:","category":"page"},{"location":"man/fitting/","page":"Fitting","title":"Fitting","text":"fitted_dens, backend_result = fit(CustomSpectrumDensity, h_sub, p0, lower_bounds, upper_bounds, shape)","category":"page"},{"location":"man/spectrum_deconvolution/#Peak-Finding-(Spectrum-Deconvolution)","page":"Peak Finding","title":"Peak Finding (Spectrum Deconvolution)","text":"","category":"section"},{"location":"man/spectrum_deconvolution/","page":"Peak Finding","title":"Peak Finding","text":"See RadiationSpectra.peakfinder","category":"page"},{"location":"man/spectrum_deconvolution/","page":"Peak Finding","title":"Peak Finding","text":"using RadiationSpectra \nh_uncal = RadiationSpectra.get_example_spectrum()\nh_decon, peakpos = RadiationSpectra.peakfinder(h_uncal)\n\nusing Plots \ngr(lw = 2.0); # hide\np_uncal = plot(h_uncal, st=:step, label=\"Uncalibrated spectrum\", c=1, lw=1.2); \np_decon = plot(peakpos, st=:vline, c=:red, label=\"Peaks\", lw=0.3);\nplot!(h_decon, st=:step, label=\"Deconvoluted spectrum\", c=1, lw=1.2); \nplot(p_uncal, p_decon, size=(800,500), layout=(2, 1)); \nsavefig(\"spectrum_devon.svg\"); # hide\nsavefig(\"spectrum_devon.pdf\"); nothing # hide","category":"page"},{"location":"man/spectrum_deconvolution/","page":"Peak Finding","title":"Peak Finding","text":"(Image: Data)","category":"page"},{"location":"LICENSE/#LICENSE","page":"LICENSE","title":"LICENSE","text":"","category":"section"},{"location":"LICENSE/","page":"LICENSE","title":"LICENSE","text":"using Markdown\nMarkdown.parse_file(joinpath(@__DIR__, \"..\", \"..\", \"LICENSE.md\"))","category":"page"},{"location":"#RadiationSpectra.jl","page":"Home","title":"RadiationSpectra.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"RadiationSpectra.jl provides tools to analyse Radiation Spectra.","category":"page"},{"location":"","page":"Home","title":"Home","text":"This includes (for now):","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"man/spectrum_deconvolution.md\"]","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"man/fitting.md\"]","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"man/spectrum_calibration.md\"]","category":"page"}]
}
