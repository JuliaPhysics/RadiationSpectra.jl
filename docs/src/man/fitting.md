# Fitting

```@contents
Pages = ["fitting.md"]
```

## Likelihood (LLH) Fit - 1D-Histogram

The underlying distribution for the Likelihood of a model to an histogram 
is the Poisson distribution. The evaluation of the likelihood of an model with a set of
parameters on the data in form of a 1D-Histogram is implemented in this package. 

Here, it is shown exemplary, how a model is defined and how it is fitted to a histgram. 

### 1. The data

The data will be in form of 1D-Histograms. 
Here, the histogram is a small subhistogram around a peak in the uncalibratum spectrum of
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


### 2. User Defined Model

Now we want to define the model, which we want to fit to the data. 

The model is a new struct which needs to be subtype of `RadiationSpectra.UvModelDensity{T}`.
Here, our model will be a Gaussian (signal) on top of flat offset (background):
```@example fitting_hist
struct MyModel{T} <: RadiationSpectra.UvModelDensity{T}
    A::T
    σ::T
    μ::T
    offset::T
end
```

Also, a method for the this type for the function `RadiationSpectra.evaluate(d::UvModelDensity, x)`
has to be defined. Note, that these are supposed to be densities as the returned value 
are multiplied to the bin width of the corresponding bin in the calculation of the likelihood. 
```@example fitting_hist
function RadiationSpectra.evaluate(d::MyModel, x)
    return d.A * pdf(Normal(d.μ, d.σ), x) + d.offset
end
```

The final step is to define a constructor for the model for a set of parameters. 
Here, the parameters are in form of s NamedTuple, which is the recommend form.
But, also a simple vector can be used. 
The package [ValueShapes.jl](https://github.com/oschulz/ValueShapes.jl) is integrated into RadiationSpectra.jl. 
Thus, also vectors or matrices can be used within NamedTuple's.
Here, only single values will be used: 

```@example fitting_hist
function MyModel(nt::NamedTuple{(:A,:μ,:σ,:offset)})
    T = promote_type(typeof.(values(nt))...)
    MyModel(T(nt.A), T(nt.σ), T(nt.μ), T(nt.offset))
end
nothing # hide
```

### 3. Set up initial guess and bounds of the parameters

Either a maximum likelihood estimation (MLE) (default) or and bayesian fit can be performed. 
In the maximum likehihood fit the package [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) is used to maximimze the likelihood. 
The bayesian fit is performed via the [BAT.jl](https://github.com/bat/BAT.jl) package.

In order to performe the MLE an initial guess, lower bounds and upper bounds for the parameter have to be given:

```@example fitting_hist
p0 = (A = 4000.0, μ = 129500, σ = 500, offset = 100) 
lower_bounds = (A = 0.0, μ = 128000, σ = 0.1, offset = 0) 
upper_bounds = (A = 30000.0, μ = 131000, σ = 1000, offset = 500) 
nothing # hide
```

Note, the package [ValueShapes.jl](https://github.com/oschulz/ValueShapes.jl) can also be used to set individual parameters constant in the fit:
```@xample
using ValueShapes
p0_fix_offset = (A = 4000.0, μ = 129500, σ = 500, offset = ConstValueShape(100)) 
```
It is enough to fix the parameter in the initial guess. 
It is not needed to also use `ConstValueShape` in the bounds. 

### 4. Performe the fit and plot the result

Now we are ready to performe the fit. 
The syntax follows the syntax of the `fit` function of 
[StatsBase.jl](https://github.com/JuliaStats/StatsBase.jl) and 
[Distributions.jl](https://github.com/JuliaStats/Distributions.jl):
`fit(::Type{Model}, data)::Model`.
Here, the Model will by our just defined model `MyModel` and the data will be the histogrom `h_sub`.
Also the initial guess and bounds have to be parsed as additional arguments:

```@example fitting_hist
fitted_model, backend_result = fit(MyModel, h_sub, p0, lower_bounds, upper_bounds)
```

The first returned argument, `fitted_model`, is an instance of `MyModel` with the fitted parameters. 

The second returned argument, `backend_result`, is, in case of an MLE fit, the returned object of the optimizer (`Optim.optimize`) used to performe the maximization. 
In case a bayesian fit, `BAT.jl` as fit backend, the second argument would be samples from BAT.jl.


### 5. Plot the fitted model:

The package provides a simple plot recipe to plot the fitted model on top of the data
through the defined `evaluate` method and the stepsize of the histogram:

```@example fitting_hist
plot(h_sub, st=:step, size=(800,400), label="Spectrum");
plot!(fitted_model, h_sub, label = "Fit");
savefig("data_hist_plus_fit.pdf"); # hide
savefig("data_hist_plus_fit.svg"); nothing # hide
```
[![Data](data_hist_plus_fit.svg)](data_hist_plus_fit.pdf)

