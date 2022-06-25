#pushfirst!(LOAD_PATH,"../src/")

using Pkg
Pkg.add(url="https://github.com/astrobc1/CurveFitParameters.jl")
Pkg.add(url="https://github.com/astrobc1/IterativeNelderMead.jl")
Pkg.add(url="https://github.com/astrobc1/EchelleBase.jl")
Pkg.add(url="https://github.com/astrobc1/EchelleReduce.jl")
Pkg.add(url="https://github.com/astrobc1/EchelleSpectralModeling.jl")
Pkg.add(url="https://github.com/astrobc1/Echelle.jl")

using Documenter
using Echelle
using EchelleBase
using EchelleReduce
using EchelleSpectralModeling

makedocs(
    sitename = "Echelle.jl",
    format = Documenter.HTML(),
    modules = [Echelle, EchelleBase, EchelleReduce, EchelleSpectralModeling],
    pages = [
        "index.md",
        "spectrographs.md",
        "reductiontutorial.md",
        "spectralmodelingtutorial.md",
        "reductionapi.md",
        "spectralmodelingapi.md",
        "modspecbehavior.md",
        "spectraldataapi.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/astrobc1/Echelle.jl.git"
)
