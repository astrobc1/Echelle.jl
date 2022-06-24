#pushfirst!(LOAD_PATH,"../src/")
pushfirst!(LOAD_PATH, "/Users/cale/Development/EchelleBase/src/")
pushfirst!(LOAD_PATH, "/Users/cale/Development/EchelleReduce/src/")
pushfirst!(LOAD_PATH, "/Users/cale/Development/CurveFitParameters/src/")
pushfirst!(LOAD_PATH, "/Users/cale/Development/IterativeNelderMead/src/")
pushfirst!(LOAD_PATH, "/Users/cale/Development/EchelleSpectralModeling/src/")
pushfirst!(LOAD_PATH, "/Users/cale/Development/Echelle/src/")

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
