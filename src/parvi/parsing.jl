using EchelleBase
using Infiltrator

function SpectralData.read_spec1d(data::SpecData1d{:parvi}, sregion::SpecRegion1d)
    f = FITSIO.FITS(data.fname)
    d = FITSIO.read(f[1])
    oi = sregion.order - ordermin(data) + 1
    λ = get_λsolution_estimate_order(sregion.order, 1)
    flux = d[oi, :, 1]
    fluxerr = d[oi, :, 2]
    mask = d[oi, :, 3]
    data.data.λ = λ
    data.data.flux = flux
    data.data.fluxerr = fluxerr
    data.data.mask = mask
    mask!(data, sregion)
    normalize!(data)
end

function SpectralData.read_image(d::MasterCal2d{:parvi}, scale_to_itime=false)
    
    image = Float64.(read(FITS(d.fname)[1]))

    # Scale slope to itime
    if scale_to_itime
        image .*= parse_itime(d)
    end
    
    return image
end

function SpectralData.read_image(d::SpecData2d{:parvi}, scale_to_itime=false)
    
    image = Float64.(read(FITS(d.fname)[1]))

    # Scale slope to itime
    if scale_to_itime
        image .*= parse_itime(d)
    end
    
    return image
end

parse_fiber_nums(d::SpecData{:parvi}) = sort!(digits(parse(Int, d.header["FIBER"])))
parse_fiber_nums(d::MasterCal2d{:parvi}) = parse_fiber_nums(d.group[1])
SpectralData.read_header(d::SpecData{:parvi}) = FITSIO.read_header(FITS(d.fname)[1])
SpectralData.parse_itime(d::SpecData{:parvi}) = d.header["EXPTIME"]
SpectralData.parse_itime(d::MasterCal2d{:parvi}) = d.group[1].header["EXPTIME"]
SpectralData.parse_object(d::SpecData{:parvi}) = d.header["OBJECT"]

function SpectralData.parse_utdate(d::SpecData{:parvi})
    t_start = parse_exposure_start_time(data)
    y, m, d = year(t_start), month(t_start), day(t_start)
    if length(m) == 1
        m = "0$m"
    end
    if length(d) == 1
        d = "0$d"
    end
    utdate = join((y, m, d))
    return utdate
end

function SpectralData.parse_sky_coord(d::SpecData{:parvi})
    if !isnothing(d.header["P200RA"]) && !isnothing(d.header["P200DEC"])
        coord = ICRSCoords(hms2rad(d.header["P200RA"]), dms2rad(d.header["P200DEC"]))
    elseif !isnothing(d.header["RA"]) && !isnothing(d.header["DEC"])
        coord = ICRSCoords(hms2rad(d.header["RA"]), dms2rad(d.header["DEC"]))
    else
        coord = ICRSCoords(0.0, 0.0)
    end
    return coord
end

SpectralData.parse_exposure_start_time(d::SpecData{:parvi}) = parse(Float64, d.header["TIMEI00"]) / 1E9 / 86400 + 2440587.5

SpectralData.get_exposure_midpoint(d::SpecData{:parvi}) = parse_exposure_start_time(d) + parse_itime(d) / (2 * 86400)
SpectralData.get_barycentric_corrections(d::SpecData{:parvi}; star_name::String) = compute_barycentric_corrections(d, star_name, "palomar")

function format_for_gfit(fname, fname_out)
    nothing
end