_your zenodo badge here_

# Submitted to JAMES 2024

**E3SM-GCAM: A synchronously coupled human coponent in the E3SM Earth system model enables novel human-Earth feedback research**

Alan V. Di Vittorio<sup>1\*</sup>, Eva Sinha<sup>2</sup>, Dalei Hao<sup>2</sup>, Balwinder Singh<sup>2</sup>, Katherine V. Calvin<sup>3</sup>,  Tim Shippert<sup>2</sup>, Pralit Patil<sup>3</sup>, and Ben Bond-Lamberty<sup>3</sup>

<sup>1 </sup> Climate and Ecosystem Sciences Division, Lawrence Berkeley National Laboratory, Berkeley, CA, United States

<sup>2 </sup> Atmospheric, Climate, & Earth Sciences Division, Pacific Northwest National Laboratory, Richland, WA, United States

<sup>3 </sup> Joint Global Change Research Institute, Pacific Northwest National Laboratory, College Park, MD, United States

\* corresponding author:  avdivittorio@lbl.gov

## Abstract

Modeling human-environment feedbacks is critical for assessing the effectiveness of climate change mitigation and adaptation strategies under a changing climate. The Energy Exascale Earth System Model (E3SM) now includes a human component, with the Global Change Analysis Model (GCAM) at its core, that is synchronously coupled with the land and atmosphere components through the E3SM coupling software. Terrestrial productivity is passed from E3SM to GCAM to make climate-responsive land use and CO2 emission projections for the next five-year period, which are interpolated and passed to E3SM annually. Key variables affected by the incorporation of these feedbacks include land use/cover change, crop prices, terrestrial carbon, local surface temperature, and climate extremes. Regional differences are more pronounced than global differences because the effects are driven primarily by differences in land use. This novel system enables a new type of scenario development and provides a modeling framework that is readily advanceable.

## Journal reference
Di Vittorio A. V., E. Sinha, D. Hao, B. Singh, K. V. Calvin, T. Shippert, P. Patel, and B. Bond-Lamberty “E3SM-GCAM: A synchronously coupled human component in the E3SM Earth system model enables novel human-Earth feedback research". JAMES, In Prep. DOI: xxxx

## Supplemental information
All supplemental materials are included in this metarepo. Paper figures and data are in `paper_figures_data`. Supplemental figures and data are in `supplemental_figures_data`. Supplemental tables are in `supplemental_tables`.


## Code reference
This is the code used for the simulations in this study.

| Model | Version | Repository Link |
|-------|---------|-----------------|
| E3SM-GCAM	| 2.1	| <https://github.com/E3SM-Project/E3SM/releases/tag/archive%2Fsinghbalwinder%2FJAMES_2024_E3SM-GCAM_Coupling_Manuscript> |

## Data reference

### Output data
The E3SM-GCAM model output data are available at: <https://portal.nersc.gov/archive/home/e/esinha/www/E3SM_GCAM_JAMES_2024>. Each future simulation has over one terrabyte of data, and the historical simulation has over two terrabytes of data.

## Contributing modeling software
These are the independent codes used in the coupled version listed above.

| Model | Version | Repository Link |
|-------|---------|-----------------|
| E3SM	| 2.1 | <https://github.com/E3SM-Project/E3SM> |
| GCAM | 6.0 | <https://github.com/JGCRI/gcam-core> |


## Reproduce the experiment
Reproducing the experiment requires access to a high performance computing cluster and explicit installation of the modeling software, which are not readily available to the average reader. The E3SM-GCAM model outputs are available as noted above, with each simulation having over one terrabyte of data. The run scripts that were used to generate model outputs are in the `run_scripts` directory. The data relevant to this  paper have been extracted from the model outputs and are provided in the `workflow` directory along with diagnostic scripts that produce the figures.


## Reproduce the figures
1. Run the following scripts in the `workflow/e3sm_analysis`, `workflow/icclim`, and `workflow/scenario_figures` directories to generate figures.

The Python scripts below have been developed using Python 3 and require installation of the modules that are imported at the beginning of each script.

The R scripts below have been developed using R 4.3.2. Required libraries are loaded at the beginning of each script using the `library()` function and must be installed as needed.

Some scripts below produce many diagnostic figures and the paper/supplemental figures have been selected from these diagnostic figures.


| Script Name | Description | How to Run |
| --- | --- | --- |
| | `workflow/e3sm_analysis` |
| `plot_spatial_diff.py` | Generate difference maps of land carbon and near surface temperature | Source the script in Python |
| `plot_ts.py` | Generate maps of land carbon, near surface temperature, and land use flux | Source the script in Python |
| | `workflow/icclim`  |
| `plot_icclim_files.py` | Generate maps of climate extreme variables | Source the script in Python |
| | `workflow/scenario_analysis` |
| `plot_e3sm_gcam_scalars.r` | Generate feedback scaling value plots | Source the script in R and run the function `plot_e3sm_gcam_scalars()` |
| `plot_gcam_crop_prices.r` | Generate crop price and associated correlation plots | Source the script in R |
| `plot_gcam_land_paper.r` | Generate scenario land allocation plot and comparison with E3SM land data plot | Source the script in R |


The map of GCAM land units showing the percent of convertible land available (Figure S1) can be recreated using the associated figure data and the map generation code in `plot_e3sm_gcam_scalars.r`. First join the appropriate land availability data to the regionXglu boundary data and then generate the map.