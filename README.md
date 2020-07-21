# SC-ION
Spatiotemporal Clustering and Inference of Omics Networks (SC-ION)

# To run
Download all R files to your computer. Run all code in the START_SCION.R file. Make sure to change the working directory to the one that contains SCION. Follow the instructions in the RShiny App to run your data.

# Known issues

- This method will automatically write over files with the same name in your working directory. If you would like to compare your results, make sure to move them from your working directory or change your directory to prevent them being written over.

- If you have previously used packages that load rlang (such as tidyverse), you may receive an error when trying to install some packages (such as dtwclust) that you need to update rlang, but it cannot be done. To fix this, you need to delete the rlang folder from your R folder, and then reinstall.

# Test data
Test data are provided in the TEST.zip folder. Download and unzip the files, then run using the following settings:

- Working directory: wherever your test files are (you can copy and paste the path from Windows/Mac file explorer)

- Target matrix: target_mat_test.csv

- Regulator matrix: reg_mat_test.csv

- Target list: target_list_test.csv

- Regulator list: reg_list_test.csv

- Clustering: You can use either Temporal or Non-Temporal and compare the results.

- Clustering matrix: cluster_test.csv

- Clustering threshold: If using temporal clustering, use 0.5. If using non-temporal clustering, use 0.25. 

- Clusters file: Leave empty

- Hub connection: No

# Version History

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
