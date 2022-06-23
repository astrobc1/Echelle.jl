using Glob
using JLD2, FileIO, FITSIO
using AstroTime
using Distributed
using NaNStatistics

using EchelleBase
using EchelleReduce

using Infiltrator

export NIRSPECReduceRecipe

struct NIRSPECReduceRecipe{E<:SpectralExtractor} <: ReduceRecipe
    data_input_path::String
    output_path::String
    sregion::SpecRegion2d
    do_dark::Bool
    do_flat::Bool
    order_height::Int
    order_spacing::Int
    base_flat_field_file::Union{String, Nothing}
    extractor::E
    extract_orders::Vector{Int}
    data::Dict{String, Any}
end

function NIRSPECReduceRecipe(;data_input_path, output_path, sregion, do_dark=false, do_flat=true, order_height, order_spacing, base_flat_field_file=nothing, extractor, extract_orders=nothing, init_data=false)
    output_path = output_path * split(data_input_path, Base.Filesystem.path_separator)[end-1] * Base.Filesystem.path_separator
    if isnothing(extract_orders)
        extract_orders = [ordermin(sregion):ordermax(sregion);]
    end
    recipe = NIRSPECReduceRecipe(data_input_path, output_path, sregion, do_dark, do_flat, order_height, order_spacing, base_flat_field_file, extractor, extract_orders, Dict{String, Any}())
    if init_data
        initialize_data(recipe)
    end
    return recipe
end

function EchelleReduce.initialize_data(recipe::NIRSPECReduceRecipe)
    
    # iSHELL science files are files that contain spc or data
    sci_files = vcat(glob("*data*.fits", recipe.data_input_path), glob("*spc*.fits", recipe.data_input_path))
    sci_files = sort(sci_files)
    recipe.data["science"] = [RawSpecData2d(sci_file, "iSHELL") for sci_file ∈ sci_files]

    # Delete bad objects
    target_names = [lowercase(parse_object(d)) for d ∈ recipe.data["science"]]
    bad = findall((target_names .== "dark") .|| (target_names .== "darks") .|| (target_names .== "flat") .|| (target_names .== "QTH"))
    deleteat!(recipe.data["science"], bad)
    
    # iSHELL flats must contain flat in the filename
    flat_files = glob("*flat*.fits", recipe.data_input_path)
    if length(flat_files) > 0
        recipe.data["flats"] = [RawSpecData2d(flat_file, "iSHELL") for flat_file ∈ flat_files]
        gas_pos = [lowercase(d.header["GASCELL"]) for d ∈ recipe.data["flats"]]
        bad = findall(gas_pos .== "in")
        deleteat!(recipe.data["flats"], bad)
        flat_groups = group_flats(recipe)
        recipe.data["master_flats"] = MasterCal2d[]
        for flat_group ∈ flat_groups
            img_nums = [ishell.parse_image_num(d) for d ∈ flat_group]
            img_start, img_end = minimum(img_nums), maximum(img_nums)
            fname = "$(recipe.output_path)calib$(Base.Filesystem.path_separator)master_flat_$(parse_utdate(flat_group[1]))_imgs$(img_start)-$(img_end).fits"
            push!(recipe.data["master_flats"], MasterCal2d(fname, flat_group))
        end
    end

    # Which to extract
    recipe.data["extract"] = recipe.data["science"]

end

function get_master_flat(data, master_flats)
    ang_seps = [separation(parse_sky_coord(mflat.group[1]), parse_sky_coord(data)) for mflat ∈ master_flats]
    ang_seps ./= 90
    time_seps = [abs(parse_exposure_start_time(mflat.group[1]) - parse_exposure_start_time(data)) for mflat ∈ master_flats]
    time_seps ./= 100
    ds = sqrt.(ang_seps.^2 + time_seps.^2)
    minds_loc = argmin(ds)
    master_flat = master_flats[minds_loc]
    return master_flat
end

function EchelleReduce.trace(recipe::NIRSPECReduceRecipe; xleft=500, xright=2048-510, n_slices=20)
    traces = Dict()
    sregions = Dict()
    for order_map ∈ recipe.data["master_flats"]
        if !isnothing(recipe.base_flat_field_file)
            offset = compute_vertical_order_drift(recipe, order_map)
            sregion_new = deepcopy(recipe.sregion)
            sregion_new.poly_bottom.coeffs[1] -= offset
            sregion_new.poly_top.coeffs[1] -= offset
            sregions[order_map] = sregion_new
        else
            sregions[order_map] = recipe.sregion
        end
        _traces = Tracing.trace(order_map, sregions[order_map], trace_pos_deg=2, min_order_spacing=recipe.order_spacing, xleft=xleft, xright=xright, n_slices=n_slices)
        for t ∈ _traces
            t["height"] -= 2
        end
        fname = "$(recipe.output_path)trace$(Base.Filesystem.path_separator)$(split(basename(order_map.fname), '.')[1])_order_map.jld"
        @save fname _traces
        traces[order_map] = _traces
    end
    return traces, sregions
end

function EchelleReduce.create_output_dirs(recipe::NIRSPECReduceRecipe)
    mkpath(recipe.output_path * "trace")
    mkpath(recipe.output_path * "spectra")
    mkpath(recipe.output_path * "calib")
end



function EchelleReduce.reduce(recipe::NIRSPECReduceRecipe)

    # Create the output directories
    create_output_dirs(recipe)

    # init data
    initialize_data(recipe)

    # Generate pre calibration images
    gen_master_calib_images(recipe)
    
    # Trace orders for appropriate images
    traces, sregions = trace(recipe)

    ny, nx = 2048, 2048

    # Fix in between orders
    for mflat ∈ recipe.data["master_flats"]
        flat_image = read_image(mflat)
        order_map_image = Tracing.gen_trace_image(traces[mflat], ny, nx, sregions[mflat])
        bad = findall(.~isfinite.(order_map_image))
        flat_image[bad] .= NaN
        FITSIO.fitswrite(mflat.fname, collect(transpose(flat_image)), header=mflat.group[1].header)
    end
    
    # Extract all desired images
    extract(recipe, traces, sregions)

end

    
function EchelleReduce.gen_master_calib_images(recipe::NIRSPECReduceRecipe)
    if recipe.do_flat
        for mflat ∈ recipe.data["master_flats"]
            image = gen_master_flat(mflat, master_bias=nothing, master_dark=nothing)
            image ./= maths.weighted_median(image, p=0.95)
            FITSIO.fitswrite(mflat.fname, collect(transpose(image)), header=mflat.group[1].header)
        end
    end
end

function EchelleReduce.extract_image(extractor::SpectralExtractor, data::SpecData2d{:ishell}, traces, sregion::SpecRegion2d, master_flat, master_dark, read_noise=0, extract_orders=nothing)
    if isnothing(extract_orders)
        extract_orders = [ordermin(sregion):ordermax(sregion);]
    end
    traces = [t for t ∈ traces if t["order"] ∈ extract_orders]
    data_image = read_image(data)
    pre_calibrate!(data_image; master_dark=master_dark, master_flat=master_flat)
    reduced_data = extract_image(extractor, data, data_image, sregion, traces, badpix_mask=nothing)
    return reduced_data
end

function EchelleReduce.extract(recipe::NIRSPECReduceRecipe, traces, sregions::Dict)
    
    pmap(1:length(recipe.data["extract"])) do i
        data = recipe.data["extract"][i]
        master_flat = recipe.do_flat ? get_master_flat(data, recipe.data["master_flats"]) : nothing
        master_dark = recipe.do_dark ? get_master_dark(data, recipe.data["master_darks"]) : nothing
        _traces = traces[master_flat]
        sregion = sregions[master_flat]
        read_noise = parse_itime(data) * detector["dark_current"] + detector["read_noise"]
        reduced_data = extract_image(recipe.extractor, data, _traces, sregion, master_flat, master_dark, read_noise, recipe.extract_orders)
        target = replace(parse_object(data), " " => "_")
        fname = "$(recipe.output_path)spectra$(Base.Filesystem.path_separator)$(splitext(basename(data.fname))[1])_$(target)_reduced.png"
        plot_extracted_spectrum(data, reduced_data, recipe.sregion, fname, [t for t ∈ _traces if t["order"] ∈ recipe.extract_orders])
        save_reduced_spectrum(recipe, data, reduced_data, sregion)
    end
end

function save_reduced_spectrum(recipe::NIRSPECReduceRecipe, data::SpecData2d{:ishell}, reduced_data, sregion::SpecRegion2d)
    n_orders = num_orders(sregion)
    reduced_data_out = fill(NaN, (n_orders, 2048, 3))
    k = 1
    for i=1:n_orders
        order = ordermin(sregion) + i - 1
        if order ∈ recipe.extract_orders && !isnothing(reduced_data[k])
            reduced_data_out[i, :, 1] .= reduced_data[k].spec1d
            reduced_data_out[i, :, 2] .= reduced_data[k].spec1derr
            reduced_data_out[i, :, 3] .= reduced_data[k].spec1dmask
            k += 1
        end
    end
    target = replace(parse_object(data), " " => "_")
    fname = "$(recipe.output_path)spectra$(Base.Filesystem.path_separator)$(splitext(basename(data.fname))[1])_$(target)_reduced.fits"
    FITSIO.fitswrite(fname, reduced_data_out, header=data.header)
end

function group_flats(recipe::NIRSPECReduceRecipe)
    flats = recipe.data["flats"]
    jds = [parse_exposure_start_time(f) for f ∈ flats]
    inds = sortperm(jds)
    flats = flats[inds]
    n_flats = length(flats)
    sep = 0.007
    prev_i = 1
    flat_groups = []
    for i=1:n_flats-1
        if jds[i+1] - jds[i] > sep
            push!(flat_groups, flats[prev_i:i])
            prev_i = i + 1
        end
    end
    push!(flat_groups, flats[prev_i:end])
    return flat_groups
end

function compute_vertical_order_drift(recipe, data)
    image1 = Float64.(read(FITS(recipe.base_flat_field_file)[1]))
    image1 = collect(transpose(image1))
    image2 = read_image(data)
    image1 .= maths.median_filter2d(image1, 5)
    image2 .= maths.median_filter2d(image2, 5)
    ny, nx = size(image1)
    lags = [-200:1:201;]
    n_lags = length(lags)
    n_slices = 100
    ccfs = fill(NaN, (n_lags, n_slices))
    yarr = [1:ny;]
    slices = Int.(floor.(collect(range(500, ny-500-1, length=n_slices))))
    for i=1:n_slices
        ii = slices[i]
        s1 = @views image1[:, ii] ./ maths.weighted_median(image1[:, ii], p=0.95)
        s2 = @views image2[:, ii] ./ maths.weighted_median(image2[:, ii], p=0.95)
        ccfs[:, i] .= maths.cross_correlate_interp(yarr, s1, yarr, s2, lags, kind="xc")
    end
    ccf = reshape(nanmedian(ccfs, dims=2), n_lags)
    lag_best = lags[maths.nanargmaximum(ccf)]
    return lag_best
end