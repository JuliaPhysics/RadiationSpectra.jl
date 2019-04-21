# Fitting

```@contents
Pages = ["fitting.md"]
```

## `Likelihood (LLH) Fit - 1D-Histogram`

1. Get a spectrum and find a peak to fit
```@example fitting_hist
using Plots, RadiationSpectra, StatsBase
myfont = Plots.font(12) # hide
pyplot(guidefont=myfont, xtickfont=myfont, ytickfont=myfont, legendfont=myfont, titlefont=myfont) # hide

h_uncal = RadiationSpectra.get_example_spectrum()
h_decon, peakpos = RadiationSpectra.peakfinder(h_uncal)

strongest_peak_bin_idx = StatsBase.binindex(h_uncal, peakpos[1])
strongest_peak_bin_width = StatsBase.binvolume(h_uncal, strongest_peak_bin_idx)
strongest_peak_bin_amplitude = h_uncal.weights[strongest_peak_bin_idx]

plot(h_uncal, st=:step, xlims=[peakpos[1] - strongest_peak_bin_width * 20, peakpos[1] + strongest_peak_bin_width * 20], size=(800,400), ylims=[0, strongest_peak_bin_amplitude * 1.1], fmt =:svg) 
```

2. Write a model function
```@example fitting_hist
function model(x, par::Vector{T}) where {T}
    scale::T = par[1]
    σ::T     = par[2]
    μ::T     = par[3]
    cp0::T   = par[4] 
    return @. scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2)) + cp0 
end
```

3. Set up the fit function [`RadiationSpectra.FitFunction{T}`](@ref).
The type, a model function, the dimensionalty of the the model and the number of parameters must be specified:
```@example fitting_hist
fitfunc = RadiationSpectra.FitFunction{Float64}( model, 1, 4); # 1 dimensional, 4 parameters 
set_fitranges!(fitfunc, ((peakpos[1] - 1000, peakpos[1] + 1000),) )
p0 = (
    A = strongest_peak_bin_amplitude * strongest_peak_bin_width * 2,
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
plot!(fitfunc, use_initial_parameters=true, lc=:green, label="Guess")
plot!(fitfunc, lc=:red, label="LLH Fit", fmt=:svg)
```

```@example fitting_hist
fitfunc # hide
```

## `LSQ Fit - 1D-Histogram`

To perfrom a LSQ Fit on a spectrum repeat the first three steps from the [Likelihood (LLH) Fit - 1D-Histogram](@ref).
Then, 

4. Performe the fit with the [`RadiationSpectra.lsqfit!`](@ref)-function and plot the result
```@example fitting_hist
RadiationSpectra.lsqfit!(fitfunc, h_uncal)

plot(h_uncal, st=:step, xlims=[peakpos[1] - strongest_peak_bin_width * 20, peakpos[1] + strongest_peak_bin_width * 20], size=(800,400), label="Spectrum", ylims=[0, strongest_peak_bin_amplitude * 1.1])
plot!(fitfunc, use_initial_parameters=true, lc=:green, label="Guess")
plot!(fitfunc, lc=:red, label="LSQ Fit", fmt=:svg)
```

```@example fitting_hist
fitfunc # hide
```

## `LSQ Fit - 1D-Data Arrays`

* Write a model function
```@example fitting_1D_data
using Plots, RadiationSpectra
myfont = Plots.font(12) # hide
pyplot(guidefont=myfont, xtickfont=myfont, ytickfont=myfont, legendfont=myfont, titlefont=myfont) # hide
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
fitfunc
```

* Performe the fit and plot the result
```@example fitting_1D_data
RadiationSpectra.lsqfit!(fitfunc, xdata, ydata, xdata_err, ydata_err) # xdata_err and ydata_err are optional

plot(xdata, ydata, xerr=xdata_err, yerr=ydata_err, st=:scatter, size=(800,400), label="Data")
plot!(fitfunc, use_initial_parameters=true, lc=:green, label="Guess")
plot!(fitfunc, lc=:red, label="LSQ Fit", fmt=:svg)
```


```@example fitting_1D_data
fitfunc # hide
```