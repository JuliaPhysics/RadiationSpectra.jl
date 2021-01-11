using Test

using Distributions
using StatsBase
using BAT

using RadiationSpectra

using BenchmarkTools

T = Float64
N = 10^6
true_pars = T[N, 0.5, 1]
h = fit(Histogram, rand(Normal(true_pars[2:3]...), N), nbins = 400)
    
p0, lower_bounds, upper_bounds, parshape, bounds = initial_parameter_guess(NormalPeakUvD{Float64}, h)

@btime NormalPeakUvD($p0)

par_vector = T[1, 0.5, 1, 0.1, 0.1]

@btime $parshape($par_vector) 
@btime RadiationSpectra._par_in_input_form($parshape($par_vector)) 
@btime NormalPeakUvD(RadiationSpectra._par_in_input_form($parshape($par_vector)))

hp = RadiationSpectra.HistLLHPrecalulations(h)

@btime RadiationSpectra.loglikelihood($hp, $NormalPeakUvD, $par_vector, $parshape)