module minerva

using FITSIO
using AstroTime
using AstroLib
using DataFrames
using Polynomials
using SpecialMatrices

using EchelleBase

SpectralData.get_spec_module(::SpecData{:minerva}) = minerva

const name = "MINERVA"
const observatory = "whipple"

const echelle_orders = [94, 122]

const detector = Dict{String, Any}(
    "gain" => 1.0,
    "read_noise" => 0.0,
    "dark_current" => 0.0,
    "nx" => 2048,
    "ny" => 2048,
    "detector_poly_bottom" => Polynomial([-46.07142857142874, 0.15125000000000033, -6.517857142857155e-05]),
    "detector_poly_top" => Polynomial([1975.4999999999993, 0.12291666666666691, -6.041666666666651e-05])
)

const lsfσ = [0.0015, 0.0022, 0.004]
const gascell_file = "iodine_gas_cell_minervanorth_nist.npz"

include("parsing.jl")

end