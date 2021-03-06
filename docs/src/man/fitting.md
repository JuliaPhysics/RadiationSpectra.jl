# Fitting

```@contents
Pages = ["fitting.md"]
```

## Likelihood (LLH) Fit - 1D-Histogram

1. Get a spectrum and find a peak to fit
```@example fitting_hist
using Plots, RadiationSpectra, StatsBase

h_uncal = RadiationSpectra.get_example_spectrum()
h_decon, peakpos = RadiationSpectra.peakfinder(h_uncal)

strongest_peak_bin_idx = StatsBase.binindex(h_uncal, peakpos[1])
strongest_peak_bin_width = StatsBase.binvolume(h_uncal, strongest_peak_bin_idx)
strongest_peak_bin_amplitude = h_uncal.weights[strongest_peak_bin_idx]

plot(h_uncal, st=:step, xlims=[peakpos[1] - strongest_peak_bin_width * 20, peakpos[1] + strongest_peak_bin_width * 20], size=(800,400), ylims=[0, strongest_peak_bin_amplitude * 1.1], fmt =:svg)
```

2. Write a model function
```@example fitting_hist
function model(x, par)
    scale = par[1]
    σ    = par[2]
    μ     = par[3]
    cp0   = par[4]
    return @. scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2)) + cp0
end
```

1. Set up the fit function [`RadiationSpectra.FitFunction{T}`](@ref).
The type, a model function, the dimensionalty of the the model and the number of parameters must be specified:
```@example fitting_hist
fitfunc = RadiationSpectra.FitFunction{Float64}( model, 1, 4); # 1 dimensional, 4 parameters
set_fitranges!(fitfunc, ((peakpos[1] - 1000, peakpos[1] + 1000),) )
p0 = (
    A = strongest_peak_bin_amplitude * 4,
    σ = strongest_peak_bin_width * 2,
    μ = peakpos[1],
    offset = 0
)
set_initial_parameters!(fitfunc, p0)
fitfunc
```

4. Performe the fit with the [`RadiationSpectra.llhfit!`](@ref)-function and plot the result
```@example fitting_hist
RadiationSpectra.llhfit!(fitfunc, h_uncal)

plot(h_uncal, st=:step, xlims=[peakpos[1] - strongest_peak_bin_width * 20, peakpos[1] + strongest_peak_bin_width * 20], size=(800,400), label="Spectrum", ylims=[0, strongest_peak_bin_amplitude * 1.1])
plot!(fitfunc, h_uncal, use_initial_parameters=true, lc=:green, label="Guess")
#plot!(fitfunc, h_uncal, lc=:red, label="LLH Fit\nΧ² = $(round(fitfunc.Chi2, digits=2))", fmt=:svg) #hide
plot!(fitfunc, h_uncal, lc=:red, label="LLH Fit", fmt=:svg)
```

```@example fitting_hist
fitfunc # hide
```

## LSQ Fit - 1D-Histogram

To perfrom a LSQ Fit on a spectrum repeat the first three steps from the [Likelihood (LLH) Fit - 1D-Histogram](@ref).
Then,

1. Perform the fit with the [`RadiationSpectra.lsqfit!`](@ref)-function and plot the result
```@example fitting_hist
RadiationSpectra.lsqfit!(fitfunc, h_uncal)

plot(h_uncal, st=:step, xlims=[peakpos[1] - strongest_peak_bin_width * 20, peakpos[1] + strongest_peak_bin_width * 20], size=(800,400), label="Spectrum", ylims=[0, strongest_peak_bin_amplitude * 1.1])
plot!(fitfunc, h_uncal, use_initial_parameters=true, lc=:green, label="Guess")
#plot!(fitfunc, h_uncal, lc=:red, label="LSQ Fit\nΧ² = $(round(fitfunc.Chi2, digits=2))", fmt=:svg) #hide
plot!(fitfunc, h_uncal, lc=:red, label="LSQ Fit", fmt=:svg)
```

```@example fitting_hist
fitfunc # hide
```
## Easy Fitting for Histograms Containing a Single Peak

There are predefined fit routines [`RadiationSpectra.fit_single_peak_histogram`](@ref), [`RadiationSpectra.fit_single_peak_histogram_refined`](@ref) that allow quick and easy fitting of standard distribution functions to histograms containing a single peak (for now supported `:Gauss`, `:Cauchy`, `:Gauss_pol1`, to be passed as keyword argument: `fit_function = ...`). Here, initial parameter guessing is attempted for you and you can obtain results in just one line.
```@example fitting_hist
subrange = (peakpos[1] - strongest_peak_bin_width * 20, peakpos[1] + strongest_peak_bin_width * 10)
plot(h_uncal, st=:step, xlims=subrange, size=(800,400), label="Spectrum", ylims=[0, strongest_peak_bin_amplitude * 1.1], lw=2)

fitfunc = fit_single_peak_histogram(h_uncal, subrange, fit_function = :Gauss_pol1 )
#plot!(fitfunc, label = "Fit Gaussian to the data", lw=3) #hide
plot!(fitfunc, label = "Fit Gaussian to the data\nΧ² = $(round(fitfunc.Chi2, digits=2))", lw=3)

fitfunc_refined = fit_single_peak_histogram_refined(h_uncal, subrange, fit_function = :Gauss_pol1, n_sig = 3 )
# plot!(fitfunc_refined, label = "Fit Gaussian to the data\nrefined parameter guessing", fmt=:svg, line = (3, :dash, :orange)) #hide
plot!(fitfunc_refined, label = "Fit Gaussian to the data\nrefined parameter guessing\nΧ² = $(round(fitfunc_refined.Chi2, digits=2))", fmt=:svg, line = (3, :dash, :orange))

```
## LSQ Fit - 1D-Data Arrays

* Write a model function
```@example fitting_1D_data
using Plots, RadiationSpectra
function model(x, par::Vector{T}) where {T}
    cp0::T   = par[1]
    cp1::T   = par[2]
    return @. cp0 + cp1 * x
end
```

* Create some random data
```@example fitting_1D_data
xdata = Float64[1, 2, 3, 6, 8, 12]
ydata = Float64[model(x, [-0.2, 0.7]) + (rand() - 0.5) for x in xdata]
xdata_err = Float64[0.5 for i in eachindex(xdata)]
ydata_err = Float64[1 for i in eachindex(xdata)]

plot(xdata, ydata, xerr=xdata_err, yerr=ydata_err, st=:scatter, size=(800,400), label="Data", fmt=:svg)
```

* Set up the fit function [`RadiationSpectra.FitFunction{T}`](@ref)
```@example fitting_1D_data
fitfunc = RadiationSpectra.FitFunction{Float64}( model, 1, 2 );
set_fitranges!(fitfunc, ((xdata[1], xdata[end]),) )
p0 = (
    offset = 0,
    linear_slope = 1
)
set_initial_parameters!(fitfunc, p0)
println(fitfunc)
```

* Performe the fit and plot the result
```@example fitting_1D_data
RadiationSpectra.lsqfit!(fitfunc, xdata, ydata, xdata_err, ydata_err) # xdata_err and ydata_err are optional

plot(xdata, ydata, xerr=xdata_err, yerr=ydata_err, st=:scatter, size=(800,400), label="Data")
plot!(fitfunc, use_initial_parameters=true, lc=:green, label="Guess")
plot!(fitfunc, lc=:red, label="LSQ Fit", fmt=:svg)
```

```@example fitting_1D_data
println(fitfunc) # hide
```

## Uncertainty Estimation

For all fit backends (`lsqfit!`, `llhfit!` and `batfit!`) one can estimate the uncertainties of the
fitted parameters
```@example fitting_hist
fitfunc.fitted_parameters
```

via
```@example fitting_hist
σs = RadiationSpectra.get_standard_deviations(fitfunc)
```

The estimation is under the assumption that the probability density function of each
parameter follows a normal distribution and that the parameters are uncorrelated.
The returned values correspond to 1σ of the normal distributions.
