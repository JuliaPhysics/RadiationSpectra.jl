# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [fixdoctests]
#
# for local builds.

using Documenter
using Pkg
ENV["PYTHON"]=""
Pkg.build("PyCall")
using Plots
pyplot(fmt=:svg)
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
    format = Documenter.HTML(   canonical = "https://JuliaHEP.github.io/RadiationSpectra.jl/stable/", 
                                prettyurls = !("local" in ARGS) 
    )
)

deploydocs(
    repo = "github.com/JuliaHEP/RadiationSpectra.jl.git",
)
