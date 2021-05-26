# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [nonstrict] [fixdoctests]
#
# for local builds.

using Documenter
using Plots; 
using RadiationSpectra

makedocs(
    sitename = "RadiationSpectra.jl",
    modules = [RadiationSpectra],
    pages=[
        "Home" => "index.md",
        "Manual" => Any[
            "Peak Finding" => "man/spectrum_deconvolution.md",
            "Fitting" => "man/fitting.md",
            "Calibration" => "man/spectrum_calibration.md",
        ],
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    format = Documenter.HTML(canonical = "https://JuliaPhysics.github.io/RadiationSpectra.jl/stable/", prettyurls = !("local" in ARGS)),
    linkcheck = ("linkcheck" in ARGS),
    strict = !("nonstrict" in ARGS),
)

deploydocs(
    repo = "github.com/JuliaPhysics/RadiationSpectra.jl.git",
    devbranch = "master",
    devurl = "master",
    forcepush = true,
    push_preview = true,
)
