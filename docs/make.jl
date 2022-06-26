pushfirst!(LOAD_PATH,"../src/")

using Documenter
println("HERE1")
using EchelleBase
println("HERE2")
using EchelleReduce
println("HERE3")
using EchelleSpectralModeling
println("HERE4")
using Echelle
println("HERE5")

makedocs(
    sitename = "Echelle.jl",
    format = Documenter.HTML(),
    println("HERE6")
    modules = [Echelle, EchelleBase, EchelleReduce, EchelleSpectralModeling],
    println("HERE7")
#    modules = [Echelle],
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
println("HERE8")
deploydocs(
    repo = "github.com/astrobc1/Echelle.jl.git"
)

println("HERE9")