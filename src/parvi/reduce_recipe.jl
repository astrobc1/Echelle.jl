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
    calib_output_path::String
    full_flat_files::Vector{String}
    fiber_flat_files::Dict{Int, Vector{String}}
    extract_fiber_flats::Bool
    dark_files::Vector{String}
    lfc_files::Union{Vector{String}, Nothing}
    extract_lfc::Bool
    badpix_mask_file::String
    do_dark::Bool
    do_flat::Bool
    extractor::E
end

function PARVIReduceRecipe(;utdate, target_input_paths, output_path, dark_files, full_flat_files, fiber_flat_files, extract_fiber_flats=true, lfc_files, extract_lfc=true, badpix_mask_file, do_dark, do_flat, extractor)
    targets = keys(target_input_paths)
    target_output_paths = Dict{String, String}(t => output_path * t * Base.Filesystem.path_separator for t ∈ keys(target_input_paths))
    calib_output_path = "$(output_path)calib_$(utdate)$(Base.Filesystem.path_separator)"
    return PARVIReduceRecipe(utdate, target_input_paths, target_output_paths, calib_output_path, full_flat_files, fiber_flat_files, extract_fiber_flats, dark_files, lfc_files, extract_lfc, badpix_mask_file, do_dark, do_flat, extractor)
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
    traces_fiber1, traces_fiber3, sregion_fiber1, sregion_fiber3 = get_traces(recipe, data)
    
    # Extract all desired images
    extract(recipe, data, traces_fiber1, traces_fiber3, sregion_fiber1, sregion_fiber3)

end

function EchelleReduce.create_output_dirs(recipe::PARVIReduceRecipe)
    for t ∈ values(recipe.target_output_paths)
        mkpath(t)
    end
    mkpath(recipe.calib_output_path)
end

EchelleReduce.get_trace_spacing(recipe::PARVIReduceRecipe, data::SpecData2d{:parvi}) = 20
EchelleReduce.get_trace_height(recipe::PARVIReduceRecipe, data::SpecData2d{:parvi}) = 12
EchelleReduce.get_extract_orders(recipe::PARVIReduceRecipe, data::SpecData2d{:parvi}, sregion::SpecRegion2d) = [SpectralRegions.ordermin(sregion):SpectralRegions.ordermax(sregion);]

function EchelleReduce.get_traces(recipe::PARVIReduceRecipe, data; xleft=400, xright=1848, n_slices=100, trace_pos_deg=2)

    # Trace fiber 1
    if 1 ∈ keys(recipe.fiber_flat_files)
        
        # Alias
        order_map = data["master_fiber_flat_fiber1"]
        
        # Some config
        sregion_fiber1 = get_specregion2d(recipe, order_map, 1)
        trace_spacing = get_trace_spacing(recipe, order_map)
        trace_height = get_trace_height(recipe, order_map)
        
        # Trace
        traces_fiber1 = Tracing.trace(order_map, sregion_fiber1, trace_pos_deg=trace_pos_deg, min_order_spacing=trace_spacing, xleft=xleft, xright=xright, n_slices=n_slices, fiber=1)
        
        # Use custom height
        for t ∈ traces_fiber1
            t["height"] = trace_height
        end

        # Save
        fname = "$(recipe.calib_output_path)$(split(basename(order_map.fname), '.')[1])_traces.jld"
        @save fname traces_fiber1
    else
        traces_fiber1 = nothing
        sregion_fiber1 = nothing
    end

    # Trace fiber 3
    if 3 ∈ keys(recipe.fiber_flat_files)

        # Alias
        order_map = data["master_fiber_flat_fiber3"]

        # Some config
        sregion_fiber3 = get_specregion2d(recipe, order_map, 3)
        trace_spacing = get_trace_spacing(recipe, order_map)
        trace_height = get_trace_height(recipe, order_map)

        # Trace
        traces_fiber3 = Tracing.trace(order_map, sregion_fiber3, trace_pos_deg=trace_pos_deg, min_order_spacing=trace_spacing, xleft=xleft, xright=xright, n_slices=n_slices, fiber=3)

        # Use custom height
        for t ∈ traces_fiber3
            t["height"] = trace_height
        end

        # Save
        fname = "$(recipe.calib_output_path)$(split(basename(order_map.fname), '.')[1])_traces.jld"
        @save fname traces_fiber3
    else
        traces_fiber3 = nothing
        sregion_fiber3 = nothing
    end

    # Return
    return traces_fiber1, traces_fiber3, sregion_fiber1, sregion_fiber3
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

function EchelleReduce.extract_image(recipe::PARVIReduceRecipe, data::SpecData2d{:parvi}, sregion_fiber1=nothing, sregion_fiber3=nothing, traces_fiber1=nothing, traces_fiber3=nothing, master_dark=nothing, master_full_flat=nothing, badpix_mask=nothing, extract_orders=nothing, fibers=nothing)

    # Load image
    data_image = read_image(data)

    # Pre calibrate
    pre_calibrate!(data_image; master_dark=master_dark, master_flat=master_full_flat)

    # Convert to pe from slope in adu
    itime = parse_itime(data)
    data_image .*= (detector["gain"] * itime)

    # Get read noise in PE
    read_noise = get_read_noise(data, detector["dark_current"], detector["read_noise"])

    # Extract orders
    if isnothing(extract_orders)
        extract_orders = get_extract_orders(recipe, data, sregion_fiber1)
    end

    # Fibers
    if isnothing(fibers)
        fibers = get_extract_fibers(recipe, data)
    end

    # Extract fiber 1
    if 1 ∈ fibers
        _traces_fiber1 = [t for t ∈ traces_fiber1 if t["order"] ∈ extract_orders]
        results_fiber1 = extract_image(recipe.extractor, data, data_image, sregion_fiber1, _traces_fiber1, badpix_mask=badpix_mask, read_noise=read_noise)
    else
        results_fiber1 = nothing
    end

    # Extract fiber 3
    if 3 ∈ fibers
        _traces_fiber3 = [t for t ∈ traces_fiber3 if t["order"] ∈ extract_orders]
        results_fiber3 = extract_image(recipe.extractor, data, data_image, sregion_fiber3, _traces_fiber3, badpix_mask=badpix_mask, read_noise=read_noise)
    else
        results_fiber3 = nothing
    end

    # Return
    return results_fiber1, results_fiber3

end

function EchelleReduce.get_specregion2d(recipe::PARVIReduceRecipe, data::SpecData2d{:parvi}, fiber=1)
    if fiber == 1
        return SpecRegion2d(
            pixmin=100, pixmax=1950,
            orderbottom=129, ordertop=85,
            poly_bottom=Polynomial([-96.96694214876034, 0.19834710743801653, -7.933884297520661e-5]) + 40,
            poly_top=Polynomial([1945.187134502924, 0.14471929824561403, -6.5906432748538e-5])
        )
    else
        return SpecRegion2d(
            pixmin=100, pixmax=1950,
            orderbottom=129, ordertop=85,
            poly_bottom=Polynomial([-96.96694214876034, 0.19834710743801653, -7.933884297520661e-5]) - 13 + 40,
            poly_top=Polynomial([1945.187134502924, 0.14471929824561403, -6.5906432748538e-5]) - 13
        )
    end
end

function get_extract_fibers(recipe::PARVIReduceRecipe, data::SpecData2d{:parvi})
    fibers = [1, 3]
    if occursin("fiber1", data.fname)
        fibers = [1]
    elseif occursin("fiber3", data.fname)
        fibers = [3]
    elseif startswith(data.fname, "LFC")
        fibers = [3]
    elseif occursin("master_lfc", data.fname)
        fibers = [1, 3]
    else
        fibers = [1]
    end
    return fibers
end

function EchelleReduce.extract(recipe::PARVIReduceRecipe, data, traces_fiber1=nothing, traces_fiber3=nothing, sregion_fiber1=nothing, sregion_fiber3=nothing)

    pmap(1:length(data["extract"])) do i
    #map(1:length(data["extract"])) do i
        
        # Alias data
        _data = data["extract"][i]

        # Alias cals
        master_dark = data["master_dark"]
        master_full_flat = data["master_full_flat"]

        # Alias badpix mask
        badpix_mask = data["badpix_mask"]

        # Which orders and fibers
        extract_orders = get_extract_orders(recipe, _data, sregion_fiber1)
        fibers = get_extract_fibers(recipe, _data)
        
        # Extract all orders and fibers
        results_fiber1, results_fiber3 = extract_image(recipe, _data, sregion_fiber1, sregion_fiber3, traces_fiber1, traces_fiber3, master_dark, master_full_flat, badpix_mask, extract_orders, fibers)

        # Plot
        plot_extracted_spectrum(recipe, _data, results_fiber1, results_fiber3, traces_fiber1, traces_fiber3, sregion_fiber1, sregion_fiber3, extract_orders)

        # Save
        save_reduced_spectrum(recipe, _data, results_fiber1, results_fiber3, sregion_fiber1, sregion_fiber3)

    end
end


function save_reduced_spectrum(recipe::PARVIReduceRecipe, data::MasterCal2d{:parvi}, results_fiber1=nothing, results_fiber3=nothing, sregion_fiber1=nothing, sregion_fiber3=nothing)
    fname = "$(recipe.calib_output_path)$(splitext(basename(data.fname))[1])_reduced.fits"
    f = FITS(fname, "w")
    header = deepcopy(data.group[1].header)
    if !isnothing(results_fiber1)
        reduced_data_out = get_extraction_result(recipe, data, sregion_fiber1, results_fiber1)
        write(f, reduced_data_out, header=header)
    else
        write(f, Float64[], header=header)
    end
    if !isnothing(results_fiber3)
        reduced_data_out = get_extraction_result(recipe, data, sregion_fiber3, results_fiber3)
        write(f, reduced_data_out)
    else
        write(f, Float64[])
    end
    close(f)
end

function save_reduced_spectrum(recipe::PARVIReduceRecipe, data::RawSpecData2d{:parvi}, results_fiber1=nothing, results_fiber3=nothing, sregion_fiber1=nothing, sregion_fiber3=nothing)
    target = parse_object(data)
    fname = "$(recipe.target_output_paths[target])$(splitext(basename(data.fname))[1])_reduced.fits"
    f = FITS(fname, "w")
    header = deepcopy(data.header)
    if !isnothing(results_fiber1)
        reduced_data_out = get_extraction_result(recipe, data, sregion_fiber1, results_fiber1)
        write(f, reduced_data_out, header=header)
    else
        write(f, Float64[], header=header)
    end
    if !isnothing(results_fiber3)
        reduced_data_out = get_extraction_result(recipe, data, sregion_fiber3, results_fiber3)
        write(f, reduced_data_out)
    else
        write(f, Float64[])
    end
    close(f)
end

function EchelleReduce.plot_extracted_spectrum(recipe::PARVIReduceRecipe, data::MasterCal2d{:parvi}, results_fiber1::Union{Vector, Nothing}=nothing, results_fiber3::Union{Vector, Nothing}=nothing, traces_fiber1::Union{Vector, Nothing}=nothing, traces_fiber3::Union{Vector, Nothing}=nothing, sregion_fiber1::Union{SpecRegion2d, Nothing}=nothing, sregion_fiber3::Union{SpecRegion2d, Nothing}=nothing, extract_orders::Union{AbstractVector, Nothing}=nothing)
    if isnothing(extract_orders)
        extract_orders = get_extract_orders(recipe, data, sregion_fiber1)
    end
    if !isnothing(results_fiber1)
        _traces_fiber1 = [t for t ∈ traces_fiber1 if t["order"] ∈ extract_orders]
        fname = "$(recipe.calib_output_path)$(splitext(basename(data.fname))[1])_fiber1_reduced.png"
        plot_extracted_spectrum(recipe, data, results_fiber1, fname, _traces_fiber1)
    end
    if !isnothing(results_fiber3)
        _traces_fiber3 = [t for t ∈ traces_fiber3 if t["order"] ∈ extract_orders]
        fname = "$(recipe.calib_output_path)$(splitext(basename(data.fname))[1])_fiber3_reduced.png"
        plot_extracted_spectrum(recipe, data, results_fiber3, fname, _traces_fiber3)
    end
end

# plot_extracted_spectrum(recipe::ReduceRecipe, data::SpecData2d, reduced_data::Vector, fname::String, traces::Vector)
function EchelleReduce.plot_extracted_spectrum(recipe::PARVIReduceRecipe, data::RawSpecData2d{:parvi}, results_fiber1::Union{Vector, Nothing}=nothing, results_fiber3::Union{Vector, Nothing}=nothing, traces_fiber1::Union{Vector, Nothing}=nothing, traces_fiber3::Union{Vector, Nothing}=nothing, sregion_fiber1::Union{SpecRegion2d, Nothing}=nothing, sregion_fiber3::Union{SpecRegion2d, Nothing}=nothing, extract_orders::Union{AbstractVector, Nothing}=nothing)
    target = parse_object(data)
    if isnothing(extract_orders)
        extract_orders = get_extract_orders(recipe, data, sregion_fiber1)
    end
    if !isnothing(results_fiber1)
        _traces_fiber1 = [t for t ∈ traces_fiber1 if t["order"] ∈ extract_orders]
        fname = "$(recipe.target_output_paths[target])$(splitext(basename(data.fname))[1])_fiber1_reduced.png"
        plot_extracted_spectrum(recipe, data, results_fiber1, fname, _traces_fiber1)
    end
    if !isnothing(results_fiber3)
        _traces_fiber3 = [t for t ∈ traces_fiber3 if t["order"] ∈ extract_orders]
        fname = "$(recipe.target_output_paths[target])$(splitext(basename(data.fname))[1])_fiber3_reduced.png"
        plot_extracted_spectrum(recipe, data, results_fiber3, fname, _traces_fiber3)
    end
end


function get_extraction_result(recipe::PARVIReduceRecipe, data::SpecData{:parvi}, sregion::SpecRegion2d, reduced_data::Vector)
    extract_orders = get_extract_orders(recipe, data, sregion)
    n_orders = abs(sregion.orderbottom - sregion.ordertop) + 1
    order_min = min(sregion.orderbottom, sregion.ordertop)
    order_max = max(sregion.orderbottom, sregion.ordertop)
    orders_all = [order_min:order_max;]
    reduced_data_out = fill(NaN, (n_orders, 2048, 3))
    for i=1:length(extract_orders)
        order = extract_orders[i]
        k = findfirst(order .== orders_all)
        if !isnothing(reduced_data[i]) && order ∈ extract_orders
            reduced_data_out[k, :, 1] .= reduced_data[i].spec1d
            reduced_data_out[k, :, 2] .= reduced_data[i].spec1derr
            reduced_data_out[k, :, 3] .= reduced_data[i].spec1dmask
        end
    end
    return reduced_data_out
end