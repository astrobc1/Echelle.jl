Reduction Tutorials
===================

#### Example 1 with iSHELL

In this example, we reduce spectra of Vega recorded with the iSHELL spectrograph in KGAS mode. Flat-field correction is performed.

```julia

# Imports
using EchelleReduce
using Echelle.ishell
using Polynomials

# Data path
data_input_path = "Vega_Test/"

# Output directory.
output_path = "reduction_results/"

# Extract all orders in KGAS mode
extract_orders = [212:230;]

# Region on the detector for iSHELL
sregion = SpecRegion2d(
    pixmin=200, pixmax=2048-200, # Left and right bounding pixels (spectral direction)
    orderbottom=212, ordertop=240, # Bottom and top order for this region on the detector (not necessary that orderbottom < ordertop).
    poly_bottom=Polynomial([-116.36685525376339, 0.20359022197025314, -5.9597390213793886e-05]), # The bottom bounding polynomial
    poly_top=Polynomial([1858.343750000002, 0.1634374999999993, -5.078125000000014e-05]) # The top bounding polynomial
)

# Standard Optimal extraction alg, for more information, see docs on extraction.
extractor = EmpiricalOptimalExtractor(;
    max_iterations=20,
    remove_background=true,
    background_smooth_width=21,
    oversample_profile=8,
    trace_pos_poly_deg=4,
    badpix_σ=5,
    extract_aperture=:auto
)

# Create the iSHELL specific recipe object, for more information, see iSHELL docs.
recipe = iSHELLReduceRecipe(
    data_input_path=data_input_path, output_path=output_path,
    sregion=sregion,
    base_flat_field_file=nothing,
    slit_height=28, order_spacing=30,
    do_flat=true, do_dark=false,
    extractor=extractor, extract_orders=extract_orders
)

# Reduce according to the recipe.
EchelleReduce.reduce(recipe)

```

The majority of the behavior of each recipe is unique for each spectrograph - tailored for their unique data products. For iSHELL, three output directories are created:

**spectra** - Contains the reduced spectra. See iSHELL docs on how these .fits files are formatted.
**calib** - Contains master calibration frames.
**trace** - Contains trace parameters derived from the master flat-fields.

