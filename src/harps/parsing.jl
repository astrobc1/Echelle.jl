using EchelleBase
using Infiltrator

function SpectralData.read_spec1d(data::SpecData1d{:harps}, sregion::SpecRegion1d)
    f = FITSIO.FITS(data.fname)
    λ = Float64.(read(f[2], "WAVE")[:] ./ 10)
    flux = Float64.(read(f[2], "FLUX")[:])
    #fluxerr = Float64.(read(f[2], "ERR")[:])
    mask = ones(length(flux))
    λ .= maths.doppler_shift_λ(λ, -1000 * data.header["ESO DRS BERV"])
    good = findall((λ .> sregion.λmin - 0.5) .&& (λ .< sregion.λmax + 0.5))
    data.data.λ = λ[good]
    data.data.flux = flux[good]
    #data.data.fluxerr = fluxerr[good]
    data.data.mask = mask[good]
    bad = findall(.~isfinite.(data.data.flux) .|| (data.data.mask .== 0) .|| (data.data.flux .<= 0))
    data.data.λ[bad] .= NaN
    data.data.flux[bad] .= NaN
    data.data.mask[bad] .= 0
    data.data.fluxerr = sqrt.(data.data.flux)
    normalize!(data, p=0.98)
end


SpectralData.read_header(d::SpecData{:harps}) = FITSIO.read_header(FITS(d.fname)[1])
SpectralData.parse_itime(d::SpecData{:harps}) = d.header["EXPTIME"]
SpectralData.parse_object(d::SpecData{:harps}) = d.header["OBJECT"]
SpectralData.parse_exposure_start_time(d::SpecData{:harps}) = d.header["MJD-OBS"]

SpectralData.get_exposure_midpoint(d::SpecData{:harps}) = d.header["ESO DRS BJD"]
SpectralData.get_barycentric_velocity(d::SpecData{:harps}) = d.header["ESO DRS BERV"] * 1E3
function SpectralData.get_barycentric_corrections(data::SpecData{:harps}; star_name)
    data.header["bjd"] = get_exposure_midpoint(data)
    data.header["bc_vel"] = get_barycentric_velocity(data)
end