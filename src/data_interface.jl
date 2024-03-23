
# Exports
export SpecData
export get_itime, get_object, get_utdate, get_sky_coord, get_exposure_start_time, get_image_num, get_timestamp
export read_header, read_header_key


# Abstract type for echelle spectra at any processing level
# S = symbol for instrument name
# L = level of processing
abstract type SpecData{S, L} end
spectrograph(::Type{<:SpecData{S, L}}) where {S, L} = S
level(::Type{<:SpecData{S, L}}) where {S, L} = L


# Interface parsing methods
function get_object end
function get_itime end
const get_exptime = get_itime
function get_utdate end
function get_sky_coord end
function get_exposure_start_time end
function get_image_num end
function get_airmass end
function get_timestamp end

# Base methods
Base.basename(data::SpecData) = basename(data.filename)
Base.show(io::IO, data::SpecData) = show(io, basename(data.filename))
Base.:(==)(data1::SpecData, data2::SpecData) = false
Base.:(==)(data1::T, data2::T) where{T<:SpecData} = data1.filename == data2.filename


# Util for reading FITS headers and keys
FITSIO.read_header(data::SpecData, hdu::Int=1) = read_header(data.filename, hdu)
function FITSIO.read_key(filename::String, key::Union{String, Int}; hdu::Int=1)
    v = FITSIO.FITS(filename) do f
        FITSIO.read_key(f[hdu], key)[1]
    end
    return v
end

FITSIO.read_key(data::SpecData; key::Union{String, Int}, hdu::Int=1) = read_key(data.filename, key; hdu)

