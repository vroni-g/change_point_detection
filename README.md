# Documentation
-----
This repository holds the code used for my Master's Thesis "Correcting for 
Multiple Testing in Change Point Detection of Global Vegetation Trends" at 
Friedrich-Schiller-University of Jena, Department of Geographic Information Science
to obtain my M.Sc. Geoinformatics in August 2021.

The repository is structured as an R package, which is essentially a copy of José
Cortés PerMuTe package (https://github.com/jcortesr/PerMuTe). He liberally 
provided the code to me.
In the folder `master_thesis` the remaining scripts and part of the input data and
results are stored. Below each script in the subfolders will be described shortly.

Not for all scripts the needed input data is available in this repository as
file size was too big. Those scripts include a note at the top. For example the 
whole LAI time series is stored on the HPC facilities of GIScience working group 
of Friedrich-Schiller-University Jena. The data was downloaded and preprocessed 
by Josè Cortés and co-workers at the Max Planck Institute for Biogeochemistry Jena. 
For further information please contact me.

Scripts starting with `FIG_[...].R` contain code to produce the figures used in 
my thesis.

## Modifications to PerMuTe
------------
The following additions have been made to José Cortés code:

* `array_to_matrix.R`: has been added by Josè Cortés to convert the 3d data matrix
to a 2d one for faster and easier computation
* `cusum_function.R`: implementation of recursive CUSUM for use in permutation correction 
(for recursive CUSUM see `strucchange` R-package)
* `get_stcs.R`: 
  + significance thresholds for cluster derivation were added for further options 
  of the `null_distribution` argument: `p-values`, `brownian_motion` for recursive 
  CUSUM and `brownian_bridge_increments` for OLS-MOSUM (both see `strucchange` R-package)
  + further changes for the minimum p-value implementation, instead of maximum statistic
  + Tippet Combining Function: for each cluster the minimum p-value/maximum statistic
  is recorded and the overall image minimum/maximum if `tippet` argument is set to `TRUE`
* `modified_cusum.R`: implementation of modified CUSUM, adapted from Lyubchich et 
al. 2020 (https://doi.org/10.1002/env.2591), with trimmed combinations and adaptive 
bootstrapping (for the original implementation by Lyubchich et al. 2020 see `funtimes` R-package);
this method was ultimately used for change point detection in the thesis
* `mosum_function.R`: implementation of OLS-MOSUM for use in permutation correction 
(for OLS-MOSUM see `strucchange` R-package)
* `perm_dist.R`: added derivation of Tippet Combining Function values as alternative to
STCS for corrected thresholding
* `threshold_data.r`: added extra thresholding option for Tippet Combining Function

## Master Thesis Folder
-------

### data_preprocessing

Data:
* `barren_land_ice_poly.rds`: spatial R object used to mask barren land and ice
* `perm_matrix.rds`: permutation matrix with 100 permutations, the first 30 
permutations were used in MCUSUM change point detection; the file was created as
computation was executed partially because of long computing times and high usage of HPC

Scripts:
* `create_mean_array_BU.R`: combines the yearly mean images of BU Gimms3g LAI data 
produced by Josè Cortés into a (non spatial) 3d array for the analysis
* `retrieve_qualityflags_BU.R`: combines the quality flags of all images of BU Gimms3g 
LAI data, for each pixel it delivers the percentage of time points with quality issues
* `FIG_qualityflags_map.R`: creates a global map of the quality issue percentage values

<a href="https://arxiv.org/pdf/2003.06222.pdf" class="uri">https://arxiv.org/pdf/2003.06222.pdf</a>

### monotonic_trend
