module Echelle

using FITSIO
using PyCall
using AstroLib, SkyCoords
using ImageFiltering
using DataInterpolations
using Polynomials
using NaNStatistics, StatsBase, Statistics


export SPEED_OF_LIGHT_MPS, COLORS_GADFLY_HEX, PATHSEP

const SPEED_OF_LIGHT_MPS = 299792458.0

const PATHSEP = Base.Filesystem.path_separator

const COLORS_GADFLY_HEX = [
    "#00BEFF", "#D4CA3A", "#FF6DAE", "#67E1B5", "#EBACFA",
    "#9E9E9E", "#F1988E", "#5DB15A", "#E28544", "#52B8AA"
]

include("utils.jl")
include("data_interface.jl")
include("skymath.jl")
include("barycenter.jl")


end
