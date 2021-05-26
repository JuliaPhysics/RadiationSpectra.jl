# Peak Finding (Spectrum Deconvolution)

See [`RadiationSpectra.peakfinder`](@ref)

```@example
using RadiationSpectra 
h_uncal = RadiationSpectra.get_example_spectrum()
h_decon, peakpos = RadiationSpectra.peakfinder(h_uncal)

using Plots 
gr(lw = 2.0); # hide
p_uncal = plot(h_uncal, st=:step, label="Uncalibrated spectrum", c=1, lw=1.2); 
p_decon = plot(peakpos, st=:vline, c=:red, label="Peaks", lw=0.3);
plot!(h_decon, st=:step, label="Deconvoluted spectrum", c=1, lw=1.2); 
plot(p_uncal, p_decon, size=(800,500), layout=(2, 1)); 
savefig("spectrum_devon.svg"); # hide
savefig("spectrum_devon.pdf"); nothing # hide
```
[![Data](spectrum_devon.svg)](spectrum_devon.pdf)