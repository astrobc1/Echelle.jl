using Glob
using JLD2, FileIO, FITSIO
using Distributed
using AstroTime
using NaNStatistics
using Infiltrator

using EchelleBase
using EchelleReduce

export PARVIReduceRecipe

const PATHSEP = Base.Filesystem.path_separator

struct PARVIReduceRecipe{E<:SpectralExtractor} <: ReduceRecipe
    utdate::String
    target_input_paths::Dict{String, String}
    target_output_paths::Dict{String, String}
    sci_exposure_fibers::Vector{Int}
    calib_output_path::String
    full_flat_files::Vector{String}
    fiber_flat_files::Dict{Int, Vector{String}}
    extract_fiber_flats::Bool
    dark_files::Vector{String}
    lfc_files::Union{Vector{String}, Nothing}
    lfc_cal_fibers::Vector{Int}
    extract_lfc::Bool
    badpix_mask_file::String
    sregion_fiber1::SpecRegion2d
    sregion_fiber3::SpecRegion2d
    do_dark::Bool
    do_flat::Bool
    extractor::E
    extract_orders::Vector{Int}
end

function PARVIReduceRecipe(;utdate, target_input_paths, output_path, sci_exposure_fibers=[1, 3], dark_files, full_flat_files, fiber_flat_files, extract_fiber_flats=true, lfc_files, lfc_cal_fibers=[1, 3], extract_lfc=true, badpix_mask_file, sregion_fiber1, sregion_fiber3, do_dark, do_flat, extractor, extract_orders=nothing)
    if isnothing(extract_orders)
        extract_orders = [ordermin(sregion_fiber1):ordermax(sregion_fiber1);]
    end
    targets = keys(target_input_paths)
    target_output_paths = Dict{String, String}(t => output_path * t * Base.Filesystem.path_separator for t ∈ keys(target_input_paths))
    calib_output_path = "$(output_path)calib_$(utdate)$(Base.Filesystem.path_separator)"
    return PARVIReduceRecipe(utdate, target_input_paths, target_output_paths, sci_exposure_fibers, calib_output_path, full_flat_files, fiber_flat_files, extract_fiber_flats, dark_files, lfc_files, lfc_cal_fibers, extract_lfc, badpix_mask_file, sregion_fiber1, sregion_fiber3, do_dark, do_flat, extractor, extract_orders)
end

function EchelleReduce.initialize_data(recipe::PARVIReduceRecipe)

    # Data
    data = Dict{String, Any}()

    # Sci files
    sci_files = String[]
    for t ∈ values(recipe.target_input_paths)
        sci_files = vcat(sci_files, glob("*$(recipe.utdate)*.fits", t))
    end

    # Create Echellograms from raw data
    data["science"] = [RawSpecData2d(f, "parvi") for f ∈ sci_files]
    data["fiber_flats_fiber1"] = [RawSpecData2d(f, "parvi") for f ∈ recipe.fiber_flat_files[1]]
    data["fiber_flats_fiber3"] = [RawSpecData2d(f, "parvi") for f ∈ recipe.fiber_flat_files[3]]
    data["darks"] = [RawSpecData2d(f, "parvi") for f ∈ recipe.dark_files]
    data["full_flats"] = [RawSpecData2d(f, "parvi") for f ∈ recipe.full_flat_files]

    # LFC files
    if !isnothing(recipe.lfc_files)
        data["lfc"] = [RawSpecData2d(f, "parvi") for f ∈ recipe.lfc_files]
    end

    # Master Darks
    if length(recipe.dark_files) > 0
        fname = "$(recipe.calib_output_path)master_dark_$(recipe.utdate).fits"
        data["master_dark"] = MasterCal2d(fname, data["darks"])
    end

    # Master Flats
    if length(recipe.full_flat_files) > 0
        fname = "$(recipe.calib_output_path)master_fullflat_$(recipe.utdate).fits"
        data["master_full_flat"] = MasterCal2d(fname, data["full_flats"])
    end
    
    # Master fiber flats
    fname = "$(recipe.calib_output_path)master_fiberflat_fiber1_$(recipe.utdate).fits"
    data["master_fiber_flat_fiber1"] = MasterCal2d(fname, data["fiber_flats_fiber1"])
    fname = "$(recipe.calib_output_path)master_fiberflat_fiber3_$(recipe.utdate).fits"
    data["master_fiber_flat_fiber3"] = MasterCal2d(fname, data["fiber_flats_fiber3"])

    # Master lfc files
    if !isnothing(recipe.lfc_files)
        fname = "$(recipe.calib_output_path)master_lfc_$(recipe.utdate).fits"
        data["master_lfc"] = MasterCal2d(fname, data["lfc"])
    end

    # Which to extract
    data["extract"] = deepcopy(data["science"])
    if recipe.extract_fiber_flats
        data["extract"] = vcat(data["extract"], [data["master_fiber_flat_fiber1"]],  [data["master_fiber_flat_fiber3"]])
    end
    if !isnothing(recipe.lfc_files) && recipe.extract_lfc
        data["extract"] = vcat(data["extract"], [data["master_lfc"]])
    end

    # Bad pixel mask (only one, load into memory)
    data["badpix_mask"] = 1.0 .- FITSIO.read(FITS(recipe.badpix_mask_file)[1])

    return data

end

function EchelleReduce.reduce(recipe::PARVIReduceRecipe)
    
    # Create the output directories
    create_output_dirs(recipe)

    # init data
    data = initialize_data(recipe)

    # Generate pre calibration images
    gen_master_calib_images(recipe, data)
    
    # Trace orders for appropriate images
    traces_fiber1, traces_fiber3 = trace(recipe, data)
    
    # Extract all desired images
    extract(recipe, data, traces_fiber1, traces_fiber3)

end

function EchelleReduce.create_output_dirs(recipe::PARVIReduceRecipe)
    for t ∈ values(recipe.target_output_paths)
        mkpath(t)
    end
    mkpath(recipe.calib_output_path)
end

function EchelleReduce.trace(recipe::PARVIReduceRecipe, data; xleft=400, xright=1848, n_slices=100)

    if 1 ∈ keys(recipe.fiber_flat_files)
        order_map = data["master_fiber_flat_fiber1"]
        traces_fiber1 = Tracing.trace(order_map, recipe.sregion_fiber1, trace_pos_deg=2, min_order_spacing=20, xleft=xleft, xright=xright, n_slices=n_slices, fiber=1)
        for t ∈ traces_fiber1
            t["height"] = 12
        end
        fname = "$(recipe.calib_output_path)$(split(basename(order_map.fname), '.')[1])_traces.jld"
        @save fname traces_fiber1
    else
        traces_fiber1 = nothing
    end

    if 3 ∈ keys(recipe.fiber_flat_files)
        order_map = data["master_fiber_flat_fiber3"]
        traces_fiber3 = Tracing.trace(order_map, recipe.sregion_fiber3, trace_pos_deg=2, min_order_spacing=20, xleft=xleft, xright=xright, n_slices=n_slices, fiber=3)
        for t ∈ traces_fiber3
            t["height"] = 12
        end
        fname = "$(recipe.calib_output_path)$(split(basename(order_map.fname), '.')[1])_traces.jld"
        @save fname traces_fiber3
    else
        traces_fiber3 = nothing
    end
    return traces_fiber1, traces_fiber3
end

function EchelleReduce.gen_master_calib_images(recipe::PARVIReduceRecipe, data)
    if recipe.do_dark
        image = gen_master_dark(data["master_dark"])
        FITSIO.fitswrite(data["master_dark"].fname, image, header=data["master_dark"].group[1].header)
    end
    if recipe.do_flat
        if recipe.do_dark
            master_dark = data["master_dark"]
        else
            master_dark = nothing
        end
        image = gen_master_flat(data["master_full_flat"], master_dark=master_dark)
        image ./= maths.weighted_median(image, p=0.95)
        FITSIO.fitswrite(data["master_full_flat"].fname, image, header=data["master_dark"].group[1].header)
    end

    image = gen_master_coadded_image(data["master_fiber_flat_fiber1"])
    FITSIO.fitswrite(data["master_fiber_flat_fiber1"].fname, image, header=data["master_fiber_flat_fiber1"].group[1].header)

    image = gen_master_coadded_image(data["master_fiber_flat_fiber3"])
    FITSIO.fitswrite(data["master_fiber_flat_fiber3"].fname, image, header=data["master_fiber_flat_fiber3"].group[1].header)

    if !isnothing(recipe.lfc_files)
        image = gen_master_coadded_image(data["master_lfc"])
        FITSIO.fitswrite(data["master_lfc"].fname, image, header=data["master_lfc"].group[1].header)
    end
end

function EchelleReduce.extract_image(extractor::SpectralExtractor, data::SpecData2d{:parvi}, sregion_fiber1::SpecRegion2d, sregion_fiber3::SpecRegion2d, traces_fiber1=nothing, traces_fiber3=nothing, master_dark=nothing, master_full_flat=nothing, badpix_mask=nothing, extract_orders=nothing, extract_fiber1=true, extract_fiber3=true)

    # Read image
    data_image = read_image(data)

    # Pre calibrate
    pre_calibrate!(data_image; master_dark=master_dark, master_flat=master_full_flat)

    # Convert to pe
    spec_mod = get_spec_module(data)
    gain = spec_mod.detector["gain"]
    itime = parse_itime(data)
    data_image .*= (gain * itime)

    # Read noise
    read_noise = itime * spec_mod.detector["dark_current"] + spec_mod.detector["read_noise"]

    # Extract fiber 1
    if extract_fiber1
        _traces_fiber1 = [t for t ∈ traces_fiber1 if t["order"] ∈ extract_orders]
        results_fiber1 = extract_image(extractor, data, data_image, sregion_fiber1, _traces_fiber1, badpix_mask=badpix_mask, read_noise=read_noise)
    else
        results_fiber1 = nothing
    end

    # Extract fiber 3
    if extract_fiber3
        _traces_fiber3 = [t for t ∈ traces_fiber3 if t["order"] ∈ extract_orders]
        results_fiber3 = extract_image(extractor, data, data_image, sregion_fiber3, _traces_fiber3, badpix_mask=badpix_mask, read_noise=read_noise)
    else
        results_fiber3 = nothing
    end

    # Return
    return results_fiber1, results_fiber3

end


function EchelleReduce.extract(recipe::PARVIReduceRecipe, data, traces_fiber1=nothing, traces_fiber3=nothing)

    pmap(1:length(data["extract"])) do i
    #map(1:length(data["extract"])) do i
        _data = data["extract"][i]
        extract_fiber1, extract_fiber3 = true, true
        master_dark = data["master_dark"]
        master_flat = data["master_full_flat"]
        badpix_mask = data["badpix_mask"]
        if _data == data["master_fiber_flat_fiber1"]
            extract_fiber1, extract_fiber3 = true, false
        end
        if _data == data["master_fiber_flat_fiber3"]
            extract_fiber1, extract_fiber3 = false, true
        end
        if _data ∈ data["science"]
            extract_fiber1, extract_fiber3 = (1 ∈ recipe.sci_exposure_fibers), (3 ∈ recipe.sci_exposure_fibers)
        end
        if "master_lfc" ∈ keys(data) && _data == data["master_lfc"]
            extract_fiber1, extract_fiber3 = (1 ∈ recipe.lfc_cal_fibers), (3 ∈ recipe.lfc_cal_fibers)
        end
        results_fiber1, results_fiber3 = extract_image(recipe.extractor, _data, recipe.sregion_fiber1, recipe.sregion_fiber3, traces_fiber1, traces_fiber3, master_dark, master_flat, badpix_mask, recipe.extract_orders, extract_fiber1, extract_fiber3)
        plot_extracted_spectrum(recipe, _data, results_fiber1, results_fiber3, traces_fiber1, traces_fiber3)
        save_reduced_spectrum(recipe, _data, results_fiber1, results_fiber3)
    end
end

function save_reduced_spectrum(recipe::PARVIReduceRecipe, data::MasterCal2d{:parvi}, results_fiber1=nothing, results_fiber3=nothing)
    fname = "$(recipe.calib_output_path)$(splitext(basename(data.fname))[1])_reduced.fits"
    f = FITS(fname, "w")
    header = deepcopy(data.group[1].header)
    if !isnothing(results_fiber1)
        reduced_data_out = get_extraction_result(recipe, results_fiber1)
        write(f, reduced_data_out, header=header)
    else
        write(f, Float64[], header=header)
    end
    if !isnothing(results_fiber3)
        reduced_data_out = get_extraction_result(recipe, results_fiber3)
        write(f, reduced_data_out)
    else
        write(f, Float64[])
    end
    close(f)
end

function save_reduced_spectrum(recipe::PARVIReduceRecipe, data::RawSpecData2d{:parvi}, results_fiber1=nothing, results_fiber3=nothing)
    target = parse_object(data)
    fname = "$(recipe.target_output_paths[target])$(splitext(basename(data.fname))[1])_reduced.fits"
    f = FITS(fname, "w")
    header = deepcopy(data.header)
    if !isnothing(results_fiber1)
        reduced_data_out = get_extraction_result(recipe, results_fiber1)
        write(f, reduced_data_out, header=header)
    else
        write(f, Float64[], header=header)
    end
    if !isnothing(results_fiber3)
        reduced_data_out = get_extraction_result(recipe, results_fiber3)
        write(f, reduced_data_out)
    else
        write(f, Float64[])
    end
    close(f)
end

function EchelleReduce.plot_extracted_spectrum(recipe::PARVIReduceRecipe, data::MasterCal2d{:parvi}, results_fiber1=nothing, results_fiber3=nothing, traces_fiber1=nothing, traces_fiber3=nothing)
    if !isnothing(results_fiber1)
        fname = "$(recipe.calib_output_path)$(splitext(basename(data.fname))[1])_fiber1_reduced.png"
        plot_extracted_spectrum(data, results_fiber1, recipe.sregion_fiber1, fname, [t for t ∈ traces_fiber1 if t["order"] ∈ recipe.extract_orders])
    end
    if !isnothing(results_fiber3)
        fname = "$(recipe.calib_output_path)$(splitext(basename(data.fname))[1])_fiber3_reduced.png"
        plot_extracted_spectrum(data, results_fiber3, recipe.sregion_fiber3, fname, [t for t ∈ traces_fiber3 if t["order"] ∈ recipe.extract_orders])
    end
end

function EchelleReduce.plot_extracted_spectrum(recipe::PARVIReduceRecipe, data::RawSpecData2d{:parvi}, results_fiber1=nothing, results_fiber3=nothing, traces_fiber1=nothing, traces_fiber3=nothing)
    target = parse_object(data)
    if !isnothing(results_fiber1)
        fname = "$(recipe.target_output_paths[target])$(splitext(basename(data.fname))[1])_fiber1_reduced.png"
        EchelleReduce.plot_extracted_spectrum(data, results_fiber1, recipe.sregion_fiber1, fname, [t for t ∈ traces_fiber1 if t["order"] ∈ recipe.extract_orders])
    end
    if !isnothing(results_fiber3)
        fname = "$(recipe.target_output_paths[target])$(splitext(basename(data.fname))[1])_fiber3_reduced.png"
        EchelleReduce.plot_extracted_spectrum(data, results_fiber3, recipe.sregion_fiber3, fname, [t for t ∈ traces_fiber3 if t["order"] ∈ recipe.extract_orders])
    end
end


function get_extraction_result(recipe::PARVIReduceRecipe, results)
    n_orders = abs(recipe.sregion_fiber1.ordertop - recipe.sregion_fiber1.orderbottom) + 1
    reduced_data_out = fill(NaN, (n_orders, 2048, 3))
    k = 1
    omin = ordermin(recipe.sregion_fiber1)
    for i=1:n_orders
        order = omin + i - 1
        if order ∈ recipe.extract_orders
            if !isnothing(results[k])
                reduced_data_out[i, :, 1] .= results[k].spec1d
                reduced_data_out[i, :, 2] .= results[k].spec1derr
                reduced_data_out[i, :, 3] .= results[k].spec1dmask
            end
            k += 1
        end
    end
    return reduced_data_out
end