module espresso

using FITSIO
using AstroTime
using AstroLib
using DataFrames

using EchelleBase

SpectralData.get_spec_module(::SpecData{:espresso}) = espresso
const name = "espresso"
const observatory = "cerro paranal"

include("parsing.jl")

const lsfσ = [0.0013, 0.0013, 0.0013]

end