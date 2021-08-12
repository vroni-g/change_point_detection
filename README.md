# Documentation

This repository holds the code used for my Master's Thesis "Correcting for 
Multiple Testing in Change Point Detection of Global Vegetation Trends" at 
Friedrich-Schiller-University of Jena, Department of Geographic Information Science
to obtain my M.Sc. Geoinformatics in August 2021.

The repository is structured as an R package, which is essentially a copy of José
Cortés PerMuTe package (https://github.com/jcortesr/PerMuTe). He liberally 
provided the code to me.
In the folder `master_thesis` the remaining scripts and part of the input data and
results are stored. Below each script in the subfolders will be described shortly.
Scripts starting with `FIG_[...].R` contain code to produce the figures used in 
my thesis. `BU` is the abbreviation for the Boston University Gimms3g LAI dataset 
used in the analysis (see Zhu et al. 2013 (https://doi.org/10.3390/rs5020927),
Chen et al. 2019 (https://doi.org/10.1038/s41893-019-0220-7)).

Not for all scripts the needed input data is available in this repository as
file size was too big. Those scripts include a note at the top. For example the 
whole LAI time series is stored on the HPC facilities of GIScience working group 
of Friedrich-Schiller-University Jena. The data was downloaded and preprocessed 
by Josè Cortés and co-workers at the Max Planck Institute for Biogeochemistry Jena. 

Thus this work is not fully reproducible, only scripts that use data provided in this 
repository can be immediately reproduced. For further information and questions please contact me.


## Modifications to PerMuTe

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

### monotonic_trend

Scripts:
*`parallel_MK_TCF_2d.R`: this script was provided by José Cortés and slightly modified,
it implements the permutation procedure in parallelized fashion with the `clustermq` R-package 
for use with SLURM and the HPC facilities in Jena
* `hpc_MK_TCF.R`: this script calls `parallel_MK_TCF_2d.R` and is the go-to script 
to define parameters of the permutation correction and ultimately send the jobs to the HPC,
the two scripts are currently written for a Mann Kendall Trend Test but can easily be 
adjusted to other methods
* `combine_MK_TCF_results.R`: combines the results of the Mann Kendall test with 
`hpc_MK_TCF.R` and derives Tippet Combining Function (TCF) values and their empirical
maximum statistic distribution for significance thresholding
* `FIG_map_MK_corrections.R`: creates a figure of three maps with the global image
of significant clusters without multiple testing correction, STCS correction and 
TCF correction
* `FIG_clustsize_MK_TCF_and_areas.R`: creates a plot of cluster size against TCF 
values and computes the areas of significant clusters for uncorrected as well as
STCS and TCF corrected results
* `FIG_qualityflags_MK_TCF.R`: creates plot of cluster size and median of quality
issues of TCF corrected significant clusters

### change_point_detection

Scripts:
* `hpc_mcusum_within_image_parallell.R`: within image parallelization of MCUSUM
testing for permuted images on HPC
* `hpc_original_parallell.R`: derives p-values as well as change point timing
estimates and selected AR orders of MCUSUM for the unpermuted original data set
using within image parallelization on HPC
* `derive_clusters.R`: for results of the two above scripts suprathreshold clusters
are derived
* `manual_quantiles_combination.R`: implements the manual combination correction 
as alternative for TCF for the MCUSUM results
* `breakpoint_types.R`: for significant cluster of original data subtrends are 
computed and different types of change points are assigned
* `FIG_map_correction_compare.R`: creates a figure of two maps with the global image
of significant clusters without multiple testing correction, STCS correction and 
manual combination correction
* `FIG_map_barplot_bp_types.R`: creates barplots for frequencies of change point 
types
* `FIG_overlay_CPD_MK.R`: creates a map and barplot of change point types for
overlay of Mann Kendall monotonic trend test and MCUSUM change point detection

### results

Data:
* `BU_LAI_MK_nperm_1000_al5.rds`: combined results of permutation procedure 
(1000 permutations) with empirical maximum statistic distributions and clusters 
of original data of Mann Kendall Trend Test with local alpha = 0.05 
* `BU_LAI_MK_nperm_1000_al10.rds`: combined results of permutation procedure 
(1000 permutations) with empirical maximum statistic distributions and clusters 
of original data of Mann Kendall Trend Test with local alpha = 0.1
* `BU_MCUSUM_cluster_original.rds`: the cluster data of original image of MCUSUM 
change point detection
* `MCUSUM_BU_orig_pvals_ar_locs.rds`: p-values, change point timing estimates and
AR order of original image of MCUSUM change point detection
* `quantiles_nperm30.rds`: cluster size based significance thresholds of manual 
combination procedure 
* `sig_mq_df.rds`: dataframe of significant MCUSUM pixel including their location,
change point timing estimate, subtrends, change point type and coincidence with Mann 
Kendall trend test
* `sig_nperm30_combined_stcsmq.rds`: matrix with positions of significant MCUSUM
pixel and vector of indices of significant clusters with manual combination


### other_analyses

Scripts:
* `compare_fpr_CPD_methods.R`: simulation to derive experimental false positive
rate of different change point detection methods
* `NA_percentage.R`: checks percentage of NA values within pixel LAI time series
* `MCUSUM_power_sims.R`: performs a simulation of MCUSUM power with artificial
breakpoint data with changing magnitude of change and differing time points of change
* `FIG_MCUSUM_power_sims.R`: creates figures of `MCUSUM_power_sims.R` results to
show power changes with changing magnitude and differing time points of change