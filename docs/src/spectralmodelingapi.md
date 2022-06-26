Spectral Modeling API
=====================

## Barycenter corrections

```@docs
compute_barycentric_corrections
```

## Ensembles - Primary containers

```@docs
SpectralRVEnsembleProblem
IterativeSpectralRVEnsembleProblem
get_init_parameters
load_templates!
```

## Spectral Models

```@docs
AbstractSpectralForwardModel
SpectralForwardModel
load_templates!
get_init_parameters
build
```

### Continuum

```@docs
PolyContinuum
SplineContinuum
estimate_continuum
```

### Wavelength Solution

```@docs
APrioriλSolution
SplineλSolution
PolyλSolution
```

#### Combs

```@docs
get_peaks
```

### Star

```@docs
AugmentedStar
augment_star!
```

### Tellurics

```@docs
TAPASTellurics
```

### Line Spread Function (LSF) / Instrument Profile (IP)

```@docs
HermiteLSF
```

### Gas Cells

```@docs
GasCell
```

## Objective Functions

```@docs
Chi2
RMS
compute_obj
```