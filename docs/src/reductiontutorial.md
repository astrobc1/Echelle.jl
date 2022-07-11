Reduction Tutorials
===================

#### Example 1 with iSHELL

Download data from [here](https://drive.google.com/file/d/1UYhsywi4XQr-jGOOslzKYvfJZ0caxDcV/view?usp=sharing).

In this example, we reduce spectra of Vega recorded with the iSHELL spectrograph in KGAS mode. A flat-field correction is performed but dark subtraction is skipped and instead a background subtraction is performed during extraction.

```julia

# Imports
using EchelleBase
using EchelleReduce
using Echelle.ishell

# Input and output path
data_input_path = "Vega_Test/"
mkpath("outputs")
output_path = "outputs/"

# Extra options through method overriding and multiple dispatch, uncomment to use.
# Which echelle orders to extract. By default, orders in KGAS are extracted (212-240).
# EchelleReduce.get_extract_orders(recipe::iSHELLReduceRecipe, data::SpecData2d{:ishell}, sregion::SpecRegion2d) = [225, 226]

# Optimal extraction
extractor = EmpiricalOptimalExtractor(;
            max_iterations=20,
            remove_background=true, background_smooth_width=21,
            oversample_profile=8, trace_pos_poly_deg=4,
            badpix_σ=5,
            extract_aperture=:auto
)

# Create the recipe class
recipe = iSHELLReduceRecipe(
         data_input_path=data_input_path, output_path=output_path,
         do_flat=true, do_dark=false,
         extractor=extractor
)

# Reduce
EchelleReduce.reduce(recipe)

```

The majority of the behavior of each recipe is unique for each spectrograph - tailored for their unique data products. For iSHELL, three output directories are created:

- **spectra** - Contains the reduced spectra. See iSHELL docs on how these .fits files are formatted.
- **calib** - Contains master calibration frames.
- **trace** - Contains trace parameters derived from the master flat-fields.

