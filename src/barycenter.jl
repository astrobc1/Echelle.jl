export get_barycentric_corrections, get_exposure_midpoint_time

function get_exposure_midpoint_time(data::SpecData, args...; kwargs...)
    return get_exposure_start_time(data, args...; kwargs...) + get_itime(data, args...; kwargs...) / (2 * 86400)
end

function get_barycentric_corrections(jdutc::Real; obs_name::String, star_name::String, zmeas::Real=0)
    barycorrpy = pyimport("barycorrpy")
    if lowercase(star_name) == "sun"
        bjd = jdutc
        bc_vel = barycorrpy.get_BC_vel(JDUTC=jdutc, starname=nothing, obsname=obs_name, leap_update=true, zmeas=0, SolSystemTarget="Sun")[1][1]
    else
        bjd = barycorrpy.utc_tdb.JDUTC_to_BJDTDB(JDUTC=jdutc, starname=star_name, obsname=obs_name, leap_update=true)[1][1]
        bc_vel = barycorrpy.get_BC_vel(JDUTC=jdutc, starname=star_name, obsname=obs_name, leap_update=true, zmeas=zmeas)[1][1]
    end
    return bjd, bc_vel
end