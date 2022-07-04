pushfirst!(LOAD_PATH,"../src/")
# pushfirst!(LOAD_PATH, "/Users/cale/Development/EchelleBase/src/")
# pushfirst!(LOAD_PATH, "/Users/cale/Development/EchelleReduce/src/")
# pushfirst!(LOAD_PATH, "/Users/cale/Development/EchelleSpectralModeling/src/")
# pushfirst!(LOAD_PATH, "/Users/cale/Development/Echelle/src/")
# pushfirst!(LOAD_PATH, "/Users/cale/Development/IterativeNelderMead/src/")
# pushfirst!(LOAD_PATH, "/Users/cale/Development/CurveFitParameters/src/")

println("Loading Documenter")
using Documenter
println("Loading EchelleBase")
using EchelleBase
println("Loading EchelleReduce")
using EchelleReduce
println("Loading EchelleSpectralModeling")
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