# SC-ION
Spatiotemporal Clustering and Inference of Omics Networks (SC-ION)

[![DOI](https://zenodo.org/badge/232898677.svg)](https://zenodo.org/badge/latestdoi/232898677)

# Citing SC-ION
If using SC-ION in your own work, please cite both this repository (using the DOI link above) and:

Clark, N.M., Nolan, T.M., Wang, P. et al. Integrated omics networks reveal the temporal signaling events of brassinosteroid response in Arabidopsis. Nat Commun 12, 5858 (2021). https://doi.org/10.1038/s41467-021-26165-3

# To run
Download all R files to your computer from Github. The files will be downloaded in a .zip folder. Extract the .zip folder. This is your directory which contains SC-ION. Make sure to change the working directory in the "START_SCION.R" file to the directory which contains all of the R code that you downloaded from Github. Then, run the START_SCION.R file to build the RShiny App. Follow the instructions in the RShiny App (or in this vignette) to run your data. 

# Known issues

- This method will automatically write over files with the same name in your working directory. If you would like to compare your results, make sure to move them from your working directory or change your directory to prevent them being written over.

- If you have previously used packages that load rlang (such as tidyverse), you may receive an error when trying to install some packages (such as dtwclust) that you need to update rlang, but it cannot be done. To fix this, you need to delete the rlang folder from your R folder, and then reinstall.

# Test data
Test data are provided in the TEST.zip folder. There are two folders inside of the TEST folder. 

The folder named Raw_Files can be used to test the script create_data_tables.R in the utilities folder. Running this script using these files should produce the files located in the Network_Files folder.

The files in the Network_Files folder can be run in SCION using the following settings. You can infer two separate networks, one for protein and one for phospho, by changing the regulator matrix and regulator list files.

- Working directory: wherever the Network_Files folder is (you can copy and paste the path from Windows/Mac file explorer)

- Target matrix: target_mat_RNA.csv

- Regulator matrix: reg_mat_protein.csv (for phospho, use reg_mat_phospho.csv)

- Target list: target_list_RNA.csv

- Regulator list: reg_list_protein.csv (for phospho, use reg_list_phospho.csv)

- Edge threshold: leave at the default value of 0.33

- Clustering: You can use either Temporal or Non-Temporal and compare the results.

- Clustering matrix: cluster_mat_protein.csv (for phospho, use cluster_mat_phospho.csv)

- Clustering threshold: If using temporal clustering, use 0.5. If using non-temporal clustering, use 2. 

- Clusters file: Leave empty

- Hub connection: No

You can use the included screenshots (shinyscreenshot_protein.png and shinyscreenshot_phospho.png) or the vignette to check your settings. This screenshot has the settings for temporal clustering.

# Version History

# Version 3.3 Hotfix - November 28, 2022

- set.seed() is now used in the main SCION function rather than in the ICA clustering function. This ensures the network inference and clustering returns the same results across multiple runs. If you would rather have different networks depending on the run, simply comment out the set.seed() line in SCION.R.

# Version 3.3 - November 7, 2022

- Converts gene IDs to R-complaint rownames. If you have experienced errors using tables with certain gene ID formats, this should now be fixed. Please be aware that the conversion to R-compliant rownames may change the format of your gene ID (for example: if your gene ID started with "0", now it will start with "X0"). Please be aware of these changes in any downstream processing of the networks.

# Version 3.2 - August 4, 2022

- A tutorial is now included in the SCION_tutorial.html file.

# Version 3.1 - June 8, 2022

- New, published test data from Zander et al, 2020, Nature Communications are now included. The previous test data have been removed. By incorporating published, citable test data, we hope to improve reproducibility of results as more features are added in the future.

- A new "utilities" folder has been added which contains scripts that may be used for pre-processing data tables and downstream analysis of the network. The create_data_tables.R script may be used to prepare tables for input into SC-ION. The example files in the TEST.zip folder may be used to understand how the script works. The NMS.R script is used to calculate a network motif-based importance score from NetMatchStar output from Cytoscape. 

# Version 3.0 - April 6, 2022

- Known compatibility: R 3.6.3, R 4.1.3

- Packages are now automatically installed and loaded using the pacman package. 

- Cluster expression is no longer automatically plotted, and the "Cluster Plots" folder is no longer automatically generated. This is to reduce the runtime of the clustering methods. Cluster expression plots can still be generated by uncommenting the plotting code in the ica_clustering.R (non-temporal) and/or dtw_clustering.R (temporal) files.

- The help text present in the RShiny App Window is now also included as a .html vignette file.

- Fixes multiple issues with network inference including adding checks for enough regulators/targets when using clustering and dealing with PTM information in regulator names. We recommend users denote PTM sites in regulator names using a "." (e.g. SOX2.S35) as this is a character not usually present in gene symbols, making it easy to separate the PTM site from the gene name if desired in downstream processing.

- Networks are no longer trimmed automatically. Rather, the user chooses an edge weight cutoff to use for trimming. For more details, see the vignette.

- Changes definition of the clustering threshold for non-temporal clustering. For more details, see the vignette.

- Using the arrows on the edge weight or clustering thresholds will increment by 0.1. NOTE: There is a known issue on MacOS where the arrows will increment by twice the desired step size (e.g. 0.2 instead of 0.1). The user can always specify an exact number by typing it into the box.

- Adds a button which will save a screenshot of the window to the working directory.

# Version 2.1 - August 23, 2021

- Same as version 2.0, but archived for publication in Nature Communications: https://doi.org/10.1038/s41467-021-26165-3

# Version 2.0 full release - August 6, 2020

- No changes from beta version

# Version 2.0 beta - July 21, 2020

- Adds another clustering algorithm (ICA) more suitable for non-temporal data

- Adds options for Temporal or Non-temporal clustering 

- Adds option to upload cluster file 

- Adds option to not cluster 

- Adds option to not connect hubs 

- Adds RS.Get.Weight.Matrix_cluster which can be used to run GENIE3 in parallel. This is not implemented in the RShiny app. If you would like to use this version, simply replace all instances of RS.Get.Weight.Matrix with RS.Get.Weight.Matrix_cluster throughout SCION.

# Version 1.1 - March 27, 2020

- Fixes edge trimming formula

- Fixes issue setting working directory

- Adds optional functions for no clustering (SCION_no_clustering) and no hub connection (SCION_no_hubs). These are not yet implemented in the RShiny app.

# Version 1.0 - Jan 9, 2020 (first stable build)
