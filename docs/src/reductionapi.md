Reduction API
=============

# Reduction Recipes

```@docs
ReduceRecipe
```

# Primary reduction methods

```@docs
EchelleReduce.reduce
initialize_data
create_output_dirs
gen_master_calib_images
get_master_flat
get_master_dark
get_traces
get_trace_spacing
get_trace_height
get_specregion2d
get_extract_orders
get_read_noise
extract
```

# Pre-extraction Calibration

```@docs
gen_master_dark
gen_master_flat
pre_calibrate!
gen_master_coadded_image
```

# Order Tracing

```@docs
trace
EchelleReduce.Tracing.gen_trace_image
```

# Spectral Extraction

```@docs
SpectralExtractor
extract_image
extract_trace
gen_model2d
plot_extracted_spectrum
```

## Spectral Extraction Utility Methods

```@docs
EchelleReduce.Extract.estimate_snr
EchelleReduce.Extract.compute_trace_positions_centroids
EchelleReduce.Extract.refine_trace_window
EchelleReduce.Extract.compute_background_1d
EchelleReduce.Extract.flag_pixels2d!
EchelleReduce.Extract.fix_bad_pixels_interp
```

## Optimal Spectral Extraction

```@docs
EmpiricalOptimalExtractor
optimal_extraction
```