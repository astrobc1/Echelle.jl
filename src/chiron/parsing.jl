using EchelleBase
using EchelleSpectralModeling
using Infiltrator

function SpectralData.read_spec1d(data::SpecData1d{:chiron}, sregion::SpecRegion1d)
    f = FITSIO.FITS(data.fname)
    d = FITSIO.read(f[1])
    oi = (echelle_orders[2] - echelle_orders[1]) - (sregion.order - echelle_orders[1]) + 1
    λ = Float64.(d[1, :, oi] ./ 10)
    flux = Float64.(d[2, :, oi])
    fluxerr = fill(1E-3, length(flux))
    mask = ones(length(flux))
    data.data.λ = λ .- gas_cell_shifts[sregion.order]
    data.data.flux = flux
    data.data.fluxerr = fluxerr
    data.data.mask = mask
    mask!(data, sregion)
    normalize!(data)
    close(f)
end


SpectralData.read_header(d::SpecData{:chiron}) = FITSIO.read_header(FITS(d.fname)[1])
SpectralData.parse_itime(d::SpecData{:chiron}) = d.header["EXPTIME"]
SpectralData.parse_object(d::SpecData{:chiron}) = d.header["OBJECT"]

SpectralData.parse_exposure_start_time(d::SpecData{:chiron}) = value(julian(Epoch(d.header["DATE-OBS"] * " TAI")))
SpectralData.get_exposure_midpoint(d::SpecData{:chiron}) = parse_exposure_start_time(d) + parse_itime(d) / (2 * 86400)
SpectralData.get_barycentric_corrections(d::SpecData{:chiron}; star_name::String) = compute_barycentric_corrections(d, star_name, "ctio")