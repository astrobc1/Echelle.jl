pushfirst!(LOAD_PATH,"../src/")

using Documenter
using EchelleBase
using EchelleReduce
using EchelleSpectralModeling
using Echelle

makedocs(
    sitename = "Echelle.jl",
    format = Documenter.HTML(),
    modules = [EchelleBase, EchelleReduce, EchelleSpectralModeling, Echelle],
    #modules = [Echelle],
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