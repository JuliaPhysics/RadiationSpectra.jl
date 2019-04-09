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
pick 858ff93 Updates make.jl for plots in documentation
s 71110b3 Adds Plots.jl & PyPlot.jl to docs/project.toml
s 29835d6 Fixes typo
s fab5d2e Updates make.jl of docs
s 5c51944 Updates make.jl format to Documenter@v0.21
s 75dc538 Add import Base.in
s 998a6fe Adds documentation
s b3c05f7 Adds more documentation
s 6c9115a Set docs `devbranch` to "dev"
s 04421ba Docs: Sets devbranch to "master" and devurl to "master"

pick 858ff93 Updates make.jl for plots in documentation
pick de78acd Adds new simple detector type `CGD` (Cuboid Germanium Detector) for cartesian coordinates. Also,o adds CartesianGrid3D method for CGD
pick 71110b3 Adds Plots.jl & PyPlot.jl to docs/project.toml
pick 29835d6 Fixes typo
pick fab5d2e Updates make.jl of docs
pick 5c51944 Updates make.jl format to Documenter@v0.21
pick d08f6fa further work on cart. coords
pick 75dc538 Add import Base.in
pick 998a6fe Adds documentation
pick 4b93f7f updates doc
pick b3c05f7 Adds more documentation
pick 6c9115a Set docs `devbranch` to "dev"
pick 04421ba Docs: Sets devbranch to "master" and devurl to "master"


pick de78acd Adds new simple detector type `CGD` (Cuboid Germanium Detector) for cartesian coordinates. Also,o adds CartesianGrid3D method for CGD
s d08f6fa further work on cart. coords
s 4b93f7f updates doc

