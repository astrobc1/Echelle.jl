module nirspec

using FITSIO
using AstroTime
using AstroLib
using DataFrames
using SkyCoords
using AstroAngles

using EchelleBase

SpectralData.get_spec_module(::SpecData{:nirspec}) = nirspec
const name = "nirspec"
const observatory = "keck"

echelle_orders = []

detector = Dict{String, Any}(
    "gain" => 1.0,
    "read_noise" => 0,
    "dark_current" => 0.0,
    "nx" => 1024,
    "ny" => 1024
)

const lsfσ = [0.03, 0.045, 0.07]

include("parsing.jl")
include("wavelength.jl")
include("reduce_recipe.jl")


end