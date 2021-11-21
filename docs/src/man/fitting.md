# Fitting

```@contents
Pages = ["fitting.md"]
```

## Likelihood (LLH) Fit - 1D-Histogram

The underlying distribution for the Likelihood of a model to a histogram 
is the Poisson distribution. The evaluation of the likelihood of an model with a set of
parameters on the data in form of a 1D-Histogram is implemented in this package. 

Here, it is shown exemplary, how a model is defined and how it is fitted to a histogram. 

### 1. The Data

The data will be in form of 1D-Histograms. 
In this example, the histogram is a small subhistogram around a peak in the uncalibrated spectrum of
a germanium detector. 

```@example fitting_hist
using Plots, RadiationSpectra, StatsBase, Distributions
gr(lw = 2.0); # hide

h_uncal = RadiationSpectra.get_example_spectrum()
h_decon, peakpos = RadiationSpectra.peakfinder(h_uncal)

strongest_peak_bin_idx = StatsBase.binindex(h_uncal, peakpos[1])
strongest_peak_bin_width = StatsBase.binvolume(h_uncal, strongest_peak_bin_idx)
strongest_peak_bin_amplitude = h_uncal.weights[strongest_peak_bin_idx]

h_sub = RadiationSpectra.subhist(h_uncal, (peakpos[1] - strongest_peak_bin_width * 20, peakpos[1] + strongest_peak_bin_width * 20))

plot(h_sub, st=:step, size=(800,400));
savefig("only_data_hist.pdf"); # hide
savefig("only_data_hist.svg"); nothing # hide
```
[![Data](only_data_hist.svg)](only_data_hist.pdf)


### 2. User Defined Spectrum Density

Now we want to define the density, which we want to fit to the data. 

It is a new `struct` which needs to be subtype of `RadiationSpectra.UvSpectrumDensity{T}`.
Here, our model will be a Gaussian (signal) on top of flat offset (background):
```@example fitting_hist
struct CustomSpectrumDensity{T} <: RadiationSpectra.UvSpectrumDensity{T}
    A::T
    σ::T
    μ::T
    offset::T
end
```

Also, a method for the this type for the function `RadiationSpectra.evaluate(d::UvSpectrumDensity, x)`
has to be defined. Note, that these are supposed to be densities as the returned value 
are multiplied to the bin volume of the corresponding bin in the calculation of the likelihood. 
```@example fitting_hist
function RadiationSpectra.evaluate(d::CustomSpectrumDensity, x)
    return d.A * pdf(Normal(d.μ, d.σ), x) + d.offset
end
```

The final step is to define a constructor for the model for a set of parameters. 
Here, the parameters are in form of `NamedTuple`, which is the recommended form.
But, also a simple vector can be used. 
The package [ValueShapes.jl](https://github.com/oschulz/ValueShapes.jl) is integrated into RadiationSpectra.jl. 
Also vectors or matrices can be used within NamedTuple's.
Here, only scalars will be used: 

```@example fitting_hist
function CustomSpectrumDensity(nt::NamedTuple{(:A,:μ,:σ,:offset)})
    T = promote_type(typeof.(values(nt))...)
    CustomSpectrumDensity(T(nt.A), T(nt.σ), T(nt.μ), T(nt.offset))
end
nothing # hide
```

### 3. Set up Initial Guess and Bounds of the Parameters

Either a maximum likelihood estimation (MLE) (default) or and bayesian fit can be performed. 
In the maximum likelihood fit the package [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) is used to maximize the likelihood. 
The bayesian fit is performed via the [BAT.jl](https://github.com/bat/BAT.jl) package.

In order to perform the MLE an initial guess, lower bounds and upper bounds for the parameter have to be given:

```@example fitting_hist
p0 = (A = 4000.0, μ = 129500, σ = 500, offset = 100) 
lower_bounds = (A = 0.0, μ = 128000, σ = 0.1, offset = 0) 
upper_bounds = (A = 30000.0, μ = 131000, σ = 1000, offset = 500) 
nothing # hide
```

### 4. Perform the Fit

Now we are ready to perform the fit. 
The syntax follows the syntax of the `fit` function of 
[StatsBase.jl](https://github.com/JuliaStats/StatsBase.jl) and 
[Distributions.jl](https://github.com/JuliaStats/Distributions.jl):
`fit(::Type{Model}, data)::Model`.
Here, the Model will by our just defined model `CustomSpectrumDensity` and the data will be the histogram `h_sub`.
Also the initial guess and bounds have to be parsed as additional arguments:

```@example fitting_hist
fitted_dens, backend_result = fit(CustomSpectrumDensity, h_sub, p0, lower_bounds, upper_bounds)
```

The first returned argument, `fitted_dens`, is an instance of `CustomSpectrumDensity` with the fitted parameters. 

The second returned argument, `backend_result`, is, in case of a MLE fit, the returned object of the optimizer 
of [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) used to perform the maximization of the likelihood. 
In case of a bayesian fit, via [`BAT.jl`](https://github.com/bat/BAT.jl) as fit backend, the second argument would be samples from BAT.jl.


### 5. Plot the fitted Density:

The package provides a simple plot recipe to plot the fitted density on top of the data
through the defined `evaluate` method and the binning of the histogram:

```@example fitting_hist
plot(h_sub, st=:step, size=(800,400), label="Spectrum");
plot!(fitted_dens, h_sub, label = "Fit");
savefig("data_hist_plus_fit.pdf"); # hide
savefig("data_hist_plus_fit.svg"); nothing # hide
```
[![Data](data_hist_plus_fit.svg)](data_hist_plus_fit.pdf)


### Optional: Fixen Parameters in the Fit

The package [ValueShapes.jl](https://github.com/oschulz/ValueShapes.jl) can be used to set individual parameters constant in the fit by defining the shape of parameters. E.g.:
```@example fitting_hist
using ValueShapes
shape = NamedTupleShape(
    A = ScalarShape{Real}(),
    μ = ScalarShape{Real}(),
    σ = ScalarShape{Real}(),
    offset = ConstValueShape{Real}(p0.offset),
)
```
The shape has to be passed to the fit function:
```@example fitting_hist
fitted_dens, backend_result = fit(CustomSpectrumDensity, h_sub, p0, lower_bounds, upper_bounds, shape)
```