using EchelleBase
using EchelleSpectralModeling

SpectralData.read_header(d::SpecData{:nirspec}) = FITSIO.read_header(FITS(d.fname)[1])
SpectralData.parse_itime(d::SpecData{:nirspec}) = d.header["ITIME"]
SpectralData.parse_object(d::SpecData{:nirspec}) = d.header["OBJECT"]
SpectralData.parse_utdate(d::SpecData{:nirspec}) = join(split(d.header["DATE_OBS"], '-'))
SpectralData.parse_sky_coord(d::SpecData{:nirspec}) = ICRSCoords(hms2rad(d.header["TCS_RA"]), dms2rad(d.header["TCS_DEC"]))
SpectralData.parse_exposure_start_time(d::SpecData{:nirspec}) = d.header["TCS_UTC"] + 2400000.5
parse_image_num(d::SpecData{:nirspec}) = parse(Int, split(split(d.fname, Base.Filesystem.path_separator)[end], '.')[5])

function SpectralData.read_image(d::MasterCal2d{:nirspec})
    image = Float64.(read(FITS(d.fname)[1]))
    image = collect(transpose(image))
    return image
end

function correct_readmath!(image, data)
    if "NDR" ∈ keys(data.header)
        image ./= float(data.header["NDR"])
    end
end

function SpectralData.read_image(d::RawSpecData2d{:nirspec})
    image = Float64.(read(FITS(d.fname)[1]))
    image = collect(transpose(image))
    correct_readmath!(image, d)
    return image
end

function SpectralData.read_spec1d(data::SpecData1d{:nirspec}, sregion)
    f = FITSIO.FITS(data.fname)
    d = FITSIO.read(f[1])
    oi = sregion.order - ordermin(data) + 1
    flux = d[oi, :, 1]
    fluxerr = d[oi, :, 2]
    mask = d[oi, :, 3]
    reverse!(flux)
    reverse!(fluxerr)
    reverse!(mask)
    data.data.flux = flux
    data.data.fluxerr = fluxerr
    data.data.mask = mask
    mask!(data, sregion)
    normalize!(data, p=0.98)
    close(f)
end

SpectralData.get_exposure_midpoint(d::SpecData{:nirspec}) = parse_exposure_start_time(d) + parse_itime(d) / (2 * 86400)
SpectralData.get_barycentric_corrections(d::SpecData{:nirspec}; star_name) = compute_barycentric_corrections(d, star_name, "irtf")