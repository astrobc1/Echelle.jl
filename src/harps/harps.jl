module harps

using FITSIO
using AstroTime
using AstroLib
using DataFrames
using Polynomials
using SpecialMatrices

using EchelleBase

SpectralData.get_spec_module(::SpecData{:harps}) = harps
const name = "harps"
const observatory = "la silla observatory"

include("parsing.jl")

const lsfσ = [0.0013, 0.0013, 0.0013]

end