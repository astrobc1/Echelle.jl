using Glob
using JLD2, FileIO, FITSIO
using AstroTime
using Distributed
using NaNStatistics

using EchelleBase
using EchelleReduce

using Infiltrator

export iSHELLReduceRecipe

struct iSHELLReduceRecipe{E<:SpectralExtractor} <: ReduceRecipe
    data_input_path::String
    output_path::String
    do_dark::Bool
    do_flat::Bool
    extractor::E
end

function iSHELLReduceRecipe(;data_input_path, output_path, do_dark=false, do_flat=true, extractor::SpectralExtractor)
    output_path = output_path * split(data_input_path, Base.Filesystem.path_separator)[end-1] * Base.Filesystem.path_separator
    recipe = iSHELLReduceRecipe(data_input_path, output_path, do_dark, do_flat, extractor)
    return recipe
end

function EchelleReduce.initialize_data(recipe::iSHELLReduceRecipe)

    # Data
    data = Dict{String, Any}()
    
    # iSHELL science files are files that contain spc or data
    sci_files = vcat(glob("*data*.fits", recipe.data_input_path), glob("*spc*.fits", recipe.data_input_path))
    sci_files = sort(sci_files)
    data["science"] = [RawSpecData2d(sci_file, "iSHELL") for sci_file ∈ sci_files]

    # Remove mislabeled cals
    target_names = [lowercase(parse_object(d)) for d ∈ data["science"]]
    bad = findall((target_names .== "dark") .|| (target_names .== "darks") .|| (target_names .== "flat") .|| (target_names .== "QTH"))
    deleteat!(data["science"], bad)
    
    # iSHELL flats must contain flat in the filename
    flat_files = glob("*flat*.fits", recipe.data_input_path)
    if length(flat_files) > 0
        data["flats"] = [RawSpecData2d(flat_file, "iSHELL") for flat_file ∈ flat_files]
        gas_pos = [lowercase(d.header["GASCELL"]) for d ∈ data["flats"]]
        bad = findall(gas_pos .== "in")
        deleteat!(data["flats"], bad)
        flat_groups = group_flats(recipe, data)
        data["master_flats"] = MasterCal2d[]
        for flat_group ∈ flat_groups
            img_nums = [ishell.parse_image_num(d) for d ∈ flat_group]
            img_start, img_end = minimum(img_nums), maximum(img_nums)
            fname = "$(recipe.output_path)calib$(Base.Filesystem.path_separator)master_flat_$(parse_utdate(flat_group[1]))_imgs$(img_start)-$(img_end).fits"
            push!(data["master_flats"], MasterCal2d(fname, flat_group))
        end
    end

    # Which to extract
    data["extract"] = data["science"]
    return data

end

function EchelleReduce.get_master_flat(recipe::iSHELLReduceRecipe, data, master_flats)
    return _get_master_flat(data, master_flats)
end

function EchelleReduce.get_master_dark(recipe::iSHELLReduceRecipe, data, master_darks)
    return _get_master_dark(data, master_darks)
end

function _get_master_flat(data, master_flats)
    ang_seps = [separation(parse_sky_coord(mflat.group[1]), parse_sky_coord(data)) for mflat ∈ master_flats]
    ang_seps ./= 90
    time_seps = [abs(parse_exposure_start_time(mflat.group[1]) - parse_exposure_start_time(data)) for mflat ∈ master_flats]
    time_seps ./= 100
    ds = sqrt.(ang_seps.^2 + time_seps.^2)
    minds_loc = argmin(ds)
    master_flat = master_flats[minds_loc]
    return master_flat
end

# Default for KGAS mode, may need adjusting!
function EchelleReduce.get_specregion2d(recipe::iSHELLReduceRecipe, data::SpecData2d{:ishell})
    sregion = SpecRegion2d(pixmin=200, pixmax=2048-200, orderbottom=212, ordertop=240,
                           poly_bottom=Polynomial([-116.36685525376339, 0.20359022197025314, -5.9597390213793886e-05]),
                           poly_top=Polynomial([1858.343750000002, 0.1634374999999993, -5.078125000000014e-05]))
    return sregion
end

EchelleReduce.get_trace_spacing(recipe::iSHELLReduceRecipe, data::SpecData2d{:ishell}) = 30
EchelleReduce.get_trace_height(recipe::iSHELLReduceRecipe, data::SpecData2d{:ishell}) = 28
EchelleReduce.get_extract_orders(recipe::iSHELLReduceRecipe, data::SpecData2d{:ishell}, sregion::SpecRegion2d) = [SpectralRegions.ordermin(sregion):SpectralRegions.ordermax(sregion);]

function EchelleReduce.get_traces(recipe::iSHELLReduceRecipe, data; xleft=500, xright=2048-510, n_slices=20)

    # Store trace params and regions for each master flat
    traces = []
    sregions = []

    # Loop over master flats
    for order_map ∈ data["master_flats"]

        # Get the sregion for this master flat
        sregion = get_specregion2d(recipe, order_map)

        # Some config
        order_spacing = get_trace_spacing(recipe, order_map)
        order_height = get_trace_height(recipe, order_map)

        # Trace the orders
        _traces = Tracing.trace(order_map, sregion, trace_pos_deg=2, min_order_spacing=order_spacing, xleft=xleft, xright=xright, n_slices=n_slices)

        # Override height results from tracing alg
        for t ∈ _traces
            t["height"] = order_height
        end

        # Save trace params to jld file
        fname = "$(recipe.output_path)trace$(Base.Filesystem.path_separator)$(split(basename(order_map.fname), '.')[1])_order_map.jld"
        @save fname _traces

        # Store results
        push!(traces, _traces)
        push!(sregions, sregion)
    end

    # Return
    return traces, sregions
end


function EchelleReduce.create_output_dirs(recipe::iSHELLReduceRecipe)
    mkpath(recipe.output_path * "trace")
    mkpath(recipe.output_path * "spectra")
    mkpath(recipe.output_path * "calib")
end

function EchelleReduce.reduce(recipe::iSHELLReduceRecipe)

    # Create the output directories
    create_output_dirs(recipe)

    # init data
    data = initialize_data(recipe)

    # Generate pre calibration images
    gen_master_calib_images(recipe, data)
    
    # Trace orders for appropriate images
    traces, sregions = get_traces(recipe, data)

    # Fix in between orders
    for (i, mflat) ∈ enumerate(data["master_flats"])
        flat_image = read_image(mflat)
        order_map_image = Tracing.gen_trace_image(traces[i], 2048, 2048, sregions[i])
        bad = findall(.~isfinite.(order_map_image))
        flat_image[bad] .= NaN
        FITSIO.fitswrite(mflat.fname, collect(transpose(flat_image)), header=mflat.group[1].header)
    end
    
    # Extract
    extract(recipe, data, traces, sregions)

end

    
function EchelleReduce.gen_master_calib_images(recipe::iSHELLReduceRecipe, data)
    if recipe.do_flat
        for mflat ∈ data["master_flats"]
            image = gen_master_flat(mflat, master_bias=nothing, master_dark=nothing)
            image ./= maths.weighted_median(image, p=0.95)
            FITSIO.fitswrite(mflat.fname, collect(transpose(image)), header=mflat.group[1].header)
        end
    end
end

function EchelleReduce.extract_image(recipe::iSHELLReduceRecipe, data::SpecData2d{:ishell}, traces::Vector, sregion::SpecRegion2d, master_flat, master_dark, read_noise=0, extract_orders=nothing)
    if isnothing(extract_orders)
        extract_orders = get_extract_orders(recipe, data, sregion)
    end
    data_image = read_image(data)
    pre_calibrate!(data_image; master_dark=master_dark, master_flat=master_flat)
    traces = [t for t ∈ traces if t["order"] ∈ extract_orders]
    reduced_data = extract_image(recipe.extractor, data, data_image, sregion, traces, badpix_mask=nothing)
    return reduced_data
end

function EchelleReduce.extract(recipe::iSHELLReduceRecipe, data, traces, sregions)

    # Extract in parallel
    map(1:length(data["extract"])) do i

        # Alias data
        _data = data["extract"][i]

        # Get cals
        master_flat = get_master_flat(recipe, _data, data["master_flats"])
        master_dark = recipe.do_dark ? get_master_dark(recipe, _data, data["master_darks"]) : nothing

        # Alias trace params and region
        k = findfirst(x -> x == master_flat, data["master_flats"])
        _traces = traces[k]
        sregion = sregions[k]

        # Get read noise in PE
        read_noise = get_read_noise(_data, detector["dark_current"], detector["read_noise"])

        # Extract full image
        reduced_data = extract_image(recipe, _data, _traces, sregion, master_flat, master_dark, read_noise)

        # Reduced filename
        target = replace(parse_object(_data), " " => "_")
        fname = "$(recipe.output_path)spectra$(Base.Filesystem.path_separator)$(splitext(basename(_data.fname))[1])_$(target)_reduced.png"

        # Plot
        extract_orders = get_extract_orders(recipe, _data, sregion)
        __traces = [t for t ∈ _traces if t["order"] ∈ extract_orders]
        plot_extracted_spectrum(recipe, _data, reduced_data, fname, __traces)

        # Save .fits file
        save_reduced_spectrum(recipe, _data, sregion, reduced_data)
    end
end

function get_extraction_result(recipe::iSHELLReduceRecipe, data::SpecData{:ishell}, sregion::SpecRegion2d, reduced_data::Vector)
    extract_orders = get_extract_orders(recipe, data, sregion)
    n_orders = abs(sregion.orderbottom - sregion.ordertop) + 1
    order_min = min(sregion.orderbottom, sregion.ordertop)
    order_max = max(sregion.orderbottom, sregion.ordertop)
    orders_all = [order_min:order_max;]
    reduced_data_out = fill(NaN, (n_orders, 2048, 3))
    for i=1:length(extract_orders)
        order = order_min + i - 1
        k = findfirst(order .== orders_all)
        if !isnothing(reduced_data[i]) && order ∈ extract_orders
            reduced_data_out[k, :, 1] .= reduced_data[i].spec1d
            reduced_data_out[k, :, 2] .= reduced_data[i].spec1derr
            reduced_data_out[k, :, 3] .= reduced_data[i].spec1dmask
        end
    end
    return reduced_data_out
end

function save_reduced_spectrum(recipe::iSHELLReduceRecipe, data::SpecData2d{:ishell}, sregion::SpecRegion2d, reduced_data::Vector)
    reduced_data_out = get_extraction_result(recipe, data, sregion, reduced_data)
    target = replace(parse_object(data), " " => "_")
    fname = "$(recipe.output_path)spectra$(Base.Filesystem.path_separator)$(splitext(basename(data.fname))[1])_$(target)_reduced.fits"
    FITSIO.fitswrite(fname, reduced_data_out, header=data.header)
end

function group_flats(recipe::iSHELLReduceRecipe, data)
    flats = data["flats"]
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