# Spectrum Calibration

The function [`RadiationSpectra.calibrate_spectrum`](@ref) return a calibrated spectrum of an uncalibrated one. 
As input it needs the uncalibrated spectrum (`StatsBase::Histogram`) and a `Vector` of known peak positions like photon lines.

The calibration is based on the assumption that the calibration function is just ``f(x) = c\cdot x``.

```@example spectrum_calibration
using Plots, RadiationSpectra 
gr(lw = 2.0); # hide

h_uncal = RadiationSpectra.get_example_spectrum()
photon_lines = [609.312, 911.204, 1120.287, 1460.830, 1764.494] # keV
h_cal, h_deconv, peakPositions, threshold, c, c_precal = RadiationSpectra.calibrate_spectrum(h_uncal, photon_lines)

p_uncal = plot(h_uncal, st=:step, label="Uncalibrated spectrum"); 
p_deconv = plot(h_deconv, st=:step, label = "Deconvoluted spectrum");
hline!([threshold], label = "threshold", lw = 1.5);
p_cal = plot(h_cal, st=:step, label="Calibrated spectrum", xlabel="E / keV"); 
vline!(photon_lines, lw=0.5, color=:red, label="Photon lines");
plot(p_uncal, p_deconv, p_cal, size=(800,700), layout=(3, 1));
savefig("calibration_of_spectrum.svg"); # hide
savefig("calibration_of_spectrum.pdf"); nothing # hide
```
[![Data](calibration_of_spectrum.svg)](calibration_of_spectrum.pdf)



Zoom into one of the peaks:
```@example spectrum_calibration
p_cal = plot(h_cal, st=:step, label="Calibrated spectrum", xlabel="E / keV", xlims=[1440, 1480], size=[800, 400]); # hide
vline!([1460.830], label="Photon line"); # hide
savefig("calibration_of_spectrum_zoom.svg"); # hide
savefig("calibration_of_spectrum_zoom.pdf"); nothing # hide
```
[![Data](calibration_of_spectrum_zoom.svg)](calibration_of_spectrum_zoom.pdf)

## Algorithm 

1. Deconvolution -> Peak finding 
2. Peak Identification - Which peak corresponds to which photon line? 
3. Fit identified peaks 
4. Fit determined peak positions vs 'true' positions (photon lines) to get the calibration constant ``c``.