module parvi

using FITSIO
using AstroTime
using AstroLib
using DataFrames
using Polynomials
using SpecialMatrices

using EchelleBase

SpectralData.get_spec_module(::SpecData{:parvi}) = parvi
const name = "PARVI"
const observatory = "palomar"

# Orders
SpectralData.orderbottom(data::SpecData{:parvi}) = 129
SpectralData.ordertop(data::SpecData{:parvi}) = 85

const detector = Dict{String, Any}(
    "gain" => 1.0,
    "read_noise" => 0.0,
    "dark_current" => 0.0,
    "nx" => 2048,
    "ny" => 2048,
    "detector_poly_bottom" => Polynomial([-46.07142857142874, 0.15125000000000033, -6.517857142857155e-05]),
    "detector_poly_top" => Polynomial([1975.4999999999993, 0.12291666666666691, -6.041666666666651e-05])
)

const lsfσ_guess_1dextract = [0.0095, 0.0095, 0.0095]

include("parsing.jl")
include("reduce_recipe.jl")
include("wavelength.jl")

end