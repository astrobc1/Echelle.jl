module chiron

using FITSIO
using AstroTime
using AstroLib
using DataFrames
using Polynomials
using SpecialMatrices

using EchelleBase

SpectralData.get_spec_module(::SpecData{:chiron}) = chiron
const name = "CHIRON"
const observatory = "ctio"

const echelle_orders = [63, 124]

# Gas cell
gascell_file = "iodine_gas_cell_chiron_master_40K.npz"

# lsf sigma
lsfσ = [0.0009, 0.0016, 0.0023]

const detector = Dict{String, Any}(
    "gain" => 1.0,
    "read_noise" => 0.0,
    "dark_current" => 0.0,
    "nx" => 2048,
    "ny" => 2048,
    "detector_poly_bottom" => Polynomial([-46.07142857142874, 0.15125000000000033, -6.517857142857155e-05]),
    "detector_poly_top" => Polynomial([1975.4999999999993, 0.12291666666666691, -6.041666666666651e-05])
)

include("parsing.jl")
include("wavelength.jl")

end