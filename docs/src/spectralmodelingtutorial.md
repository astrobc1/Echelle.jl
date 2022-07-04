Spectral Modeling Tutorials
===========================

There is currently one primary container for spectral modeling: `IterativeSpectralRVEnsembleProblem`. Regarldess of the spectrograph, the output structure is identical. A new folder is created within the provided `output_path` variable at `output_path/[spectrograph]_[tag]/`. Within this folder will be a subfolder for each order (or chunk). Within a given order/chunk folder, one finds Fits, RVs, and Templates, as well as the primary ensemble pbject stored in a `.jld` file.

The structure is summed up as follows for chunks defined by orders, where `label[M]` = `Order[M]` or `Chunk[M]`.

- `output_path/[spectrograph]_[tag]/`
- `label[M]`
    - `ensemble_Order[M].jld`
    - `Fits/`
        - `[filename]_label[M]_iter1.png`
        - `[filename]_label[M]_iter2.png`
        - `...`
        - `optimization_results_label[M].jld`
    - `RVs/`
        - `rvs_label[M]_iter1.png`
        - `rvs_label[M]_iter2.png`
        - `...`
        - `rvs_label[M].jld`
    - `Templates/`
        - `stellar_templates_label[M].txt`

# Example 1: iSHELL, KGAS mode, $^{13}CH_{4}$ gas cell

In this example, we fit two nights of spectra of Barnard's Star (GJ 699) recorded with the iSHELL spectrograph in KGAS mode.

```julia

# Imports
using EchelleBase
using EchelleSpectralModeling
using Echelle.ishell

# Config
spectrograph = "iSHELL"
data_input_path = "GJ_699_iSHELL_Example/"
filelist = "filelist_example.txt"
output_path = "/Users/cale/Research/specmodel_results/"
star_name = "GJ_699"
tag = "gj699_example"
do_orders = [219, 222, 226]
templates_path = "SpectralTemplates/"

# iSHELL gas cell depth (= 0.97, fixed).
# The deviation from unity is likely due to the off-axis angle the gas cell is placed in iSHELL vs. the FTS scan.
τ_gascell = [ishell.τ_gascell, ishell.τ_gascell, ishell.τ_gascell]

# Loop over orders
# All results will be saved to specmodel_results/ishell_gj699_example/Chunki/.../
for order ∈ do_orders

    # Pixel bounds for this order. iSHELL spectra are 2048 pixels long.
    pixmin, pixmax = 200, 1848

    # Get estimated wavelength grid in order to set corresponding wavelength bounds for this order
    λ_estimate = ishell.get_λsolution_estimate_order(order)

    # Define the region
    sregion = SpecRegion1d(pixmin=pixmin, pixmax=pixmax, λmin=λ_estimate[pixmin], λmax=λ_estimate[pixmax], order=order)

    # Create the model
    model = SpectralForwardModel(λsolution=PolyλSolution(deg=2, bounds=[-0.05, 0.05]),
                                 continuum=SplineContinuum(n_splines=6, bounds=[0.7, 1.3]),
                                 lsf=HermiteLSF(deg=0, σ_guess=ishell.lsfσ_guess_kgas_0375),
                                 star=AugmentedStar(star_name=star_name),
                                 gascell=GasCell(input_file=templates_path * ishell.gascell_file, depth_guess=τ_gascell),
                                 tellurics=TAPASTellurics(input_file=templates_path * "TAPAS_tellurics_maunakea.npz", vel_guess=[-200, 20, 200]),
                                 sregion=sregion, oversample=4)

    # Simple RMS objective (uniform weights)
    obj = RMS()

    # Simple Augmenter (weighted median)
    augmenter = WeightedMedianAugmenter()

    # Create the ensemble
    ensemble = IterativeSpectralRVEnsembleProblem(spectrograph=spectrograph, filenames, model=model, obj=obj, augmenter=augmenter)
    
    # Run RVs for this order
    compute_rvs(ensemble, output_path=output_path, tag=tag, n_iterations=5, do_ccf=false, initial_star_from_data=false, verbose=true)

end

```

Now we look at the results of the spectral forward modeling and generate RVs which are co-added across orders and/or multiple exposures iSHELL results.

```julia

```

# Example 2A: HARPS, 51-Pegasi in parallel (single computing node).

In this example, we utilize data recorded with HARPS-S downloaded from the ESO archive, which can be downloaded [here](https://drive.google.com/file/d/1-vrL2qzsWv8x_StDqOwQXJ9eqxAGFWsM/view?usp=sharing). Extract the zip once downloaded. Note this will take some time! One can limit the wavelength range to something smaller (and change the number of chunks accordingly) in the line `chunks = range(480, 700, length=101)`.

```julia

# Imports
using Distributed
using EchelleBase
using EchelleSpectralModeling
using Echelle.harps

# Parallel setup
n_cores = 2 # Increase/decrease as desired
addprocs(n_cores)
@everywhere begin
    using EchelleBase
    using EchelleSpectralModeling
    using Echelle.harps
end

# Basic info
spectrograph = "HARPS"
data_input_path = "51Peg_HARPS_ARCHIVE_EXAMPLE/"
filelist = "filelist.txt"
output_path = "specmodel_results/"
star_name = "51_Peg"
tag = "51peg_example"
templates_path = "SpectralTemplates/"

# Archival HARPS data is already stiched together, so we define our own chunks: 100 chunks from 480 - 700 nm.
chunks = range(480, 700, length=101)

# Loop over chunks
# All results will be saved to specmodel_results/harps_51peg_example/Chunki/.../
for i=1:length(chunks)-1
    
    # Define the region for this chunk. Note in this case the bounds are set by the wavelength values in the lab frame, not by pixels.
    sregion = SpecRegion1d(pixmin=nothing, pixmax=nothing, λmin=chunks[i], λmax=chunks[i+1], label="Chunk$i")

    # Create the model
    # Note we use a static lsf width, which while may not be accurate, still yields accurate relative RVs from current tests.
    # We also steal the telluric template from CTIO instead of La Silla, which is reasonable for this example as tellurics are sparse in the visible.
    model = SpectralForwardModel(λsolution=APrioriλSolution(),
                                 continuum=PolyContinuum(deg=1, coeffs_guess=Dict(0=>[0.95, 1.0, 1.1], 1=>[-1E-4, 1E-5, 1E-4])),
                                 lsf=HermiteLSF(deg=0, σ_guess=harps.lsfσ),
                                 star=AugmentedStar(star_name=star_name),
                                 tellurics=TAPASTellurics(input_file=templates_path * "TAPAS_tellurics_ctio.npz", vel_guess=[-300, 10, 300], min_feature_depth=0.1),
                                 sregion=sregion, oversample=8)

    # Simple RMS objective (uniform weights)
    obj = RMS()

    # Simple Augmenter (weighted median)
    augmenter = WeightedMedianAugmenter()

    # Create the ensemble
    ensemble = IterativeSpectralRVEnsembleProblem(spectrograph=spectrograph, filenames, model=model, obj=obj, augmenter=augmenter)
    
    # Run RVs for this order
    compute_rvs(ensemble, output_path=output_path, tag=tag, n_iterations=5, do_ccf=false, initial_star_from_data=true, verbose=true)

end

# Remove workers
rmprocs(workers())

```

Now we look at the results of the spectral forward modeling and generate RVs which are co-added across orders and/or multiple exposures from the HARPS results.

```julia



```
