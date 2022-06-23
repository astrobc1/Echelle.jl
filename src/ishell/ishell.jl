module ishell

using FITSIO
using AstroTime
using AstroLib
using DataFrames
using SkyCoords
using AstroAngles

using EchelleBase

SpectralData.get_spec_module(::SpecData{:ishell}) = ishell
const name = "iSHELL"
const observatory = "irtf"

# Orders
SpectralData.orderbottom(data::SpecData{:ishell}) = 212
SpectralData.ordertop(data::SpecData{:ishell}) = 229

const detector = Dict(
    "gain" => 1.8,
    "read_noise" => 8.0,
    "dark_current" => 0.05,
    "nx" => 2048,
    "ny" => 2048
)

# Some info for spectral forward modeling
const gascell_file = "methane_gas_cell_ishell_kgas.npz"
const gascell_depth = 0.97
const lsfσ_guess_kgas_0375 = [0.01, 0.013, 0.016]

# Parsing
include("parsing.jl")

# Default wavelength info
include("wavelength.jl")

# Reduce recipe
include("reduce_recipe.jl")


end