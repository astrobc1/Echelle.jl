using EchelleBase
using EchelleSpectralModeling

function SpectralData.read_spec1d(data::SpecData1d{:minerva}, sregion::SpecRegion1d)
    f = FITSIO.FITS(data.fname)
    d = FITSIO.read(f[1])
    oi = (echelle_orders[2] - echelle_orders[1]) - (sregion.order - echelle_orders[1]) + 1
    λ = d[1, :, oi] ./ 10
    tel = parse_telescope(data)
    if tel == 4
        λ .+= 0.0133
    end
    flux = d[2, :, oi]
    fluxerr = d[3, :, oi]
    mask = d[4, :, oi]
    data.data.λ = λ
    data.data.flux = flux
    data.data.fluxerr = fluxerr
    data.data.mask = mask
    mask!(data, sregion)
    normalize!(data, p=0.98)
end


parse_telescope(d::SpecData{:minerva}) = parse(Int, basename(d.fname)[end-5])
parse_fiber_nums(d::SpecData{:minerva}) = [parse_telescope(d)]
SpectralData.read_header(d::SpecData{:minerva}) = FITSIO.read_header(FITS(d.fname)[1])
SpectralData.parse_itime(d::SpecData{:minerva}) = d.header["EXPTIME"]
SpectralData.parse_exposure_start_time(d::SpecData{:minerva}) = d.header["JD"]

function SpectralData.get_exposure_midpoint(data)
    return parse_exposure_start_time(data) + parse_itime(data) / (2 * 86400)
end

SpectralData.get_barycentric_corrections(d::SpecData{:minerva}; star_name) = compute_barycentric_corrections(d, star_name, "whipple")

function get_λsolution_estimate_from_file(fname, order)
    oi = (echelle_orders[2] - echelle_orders[1]) - (order - echelle_orders[1]) + 1
    f = FITS(fname)
    data = read(f[1])
    λ = data[1, :, oi] ./ 10
    return λ
end

function SpectralData.get_λsolution_estimate(data::SpecData1d{:minerva}, sregion::SpecRegion1d)
    return data.data.λ
end
