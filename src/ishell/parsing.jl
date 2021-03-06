using EchelleBase
using EchelleSpectralModeling

# Header parsing
SpectralData.read_header(d::SpecData{:ishell}) = FITSIO.read_header(FITS(d.fname)[1])
SpectralData.parse_itime(d::SpecData{:ishell}) = d.header["ITIME"]
SpectralData.parse_object(d::SpecData{:ishell}) = d.header["OBJECT"]
SpectralData.parse_utdate(d::SpecData{:ishell}) = join(split(d.header["DATE_OBS"], '-'))
SpectralData.parse_sky_coord(d::SpecData{:ishell}) = ICRSCoords(hms2rad(d.header["TCS_RA"]), dms2rad(d.header["TCS_DEC"]))
SpectralData.parse_exposure_start_time(d::SpecData{:ishell}) = d.header["TCS_UTC"] + 2400000.5

# Filename parsing
parse_image_num(d::SpecData{:ishell}) = parse(Int, split(split(d.fname, Base.Filesystem.path_separator)[end], '.')[5])

# Read in master cal
function SpectralData.read_image(d::MasterCal2d{:ishell})
    image = Float64.(read(FITS(d.fname)[1]))
    image = collect(transpose(image))
    return image
end

# Correct NDRs
function correct_readmath!(image, data)
    if "NDR" ∈ keys(data.header)
        image ./= float(data.header["NDR"])
    end
end

# Read in Raw frame
function SpectralData.read_image(d::RawSpecData2d{:ishell})
    image = Float64.(read(FITS(d.fname)[1]))
    image = collect(transpose(image))
    correct_readmath!(image, d)
    return image
end

# Read in reduced spectrum
function SpectralData.read_spec1d(data::SpecData1d{:ishell}, sregion)
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
end

# Barycenter info
SpectralData.get_exposure_midpoint(d::SpecData{:ishell}) = parse_exposure_start_time(d) + parse_itime(d) / (2 * 86400)
SpectralData.get_barycentric_corrections(d::SpecData{:ishell}; star_name) = compute_barycentric_corrections(d, star_name, "irtf")