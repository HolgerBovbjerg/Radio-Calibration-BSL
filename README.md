# Radio-Calibration-BSL
This repository contains work done as part of a semester project at the Signal Procesing master at Aalborg University.  
This includes code implementations and some examples. Some scripts have included capability of parallel computation.
To run the examples the repository including subfolders must be added to the MATLAB path.
## Folder structure 
In [`/ABC`](https://github.com/HolgerBovbjerg/Radio-Calibration-BSL/tree/main/ABC) code involving the ABC method found this includes:
* Simple rejection ABC implementation - 'ABC_REJ.m'.
* PMC-ABC method proposed by Ayush Bharti - 'ABC_REJ_PMC.m'.

In [`/BSL`](https://github.com/HolgerBovbjerg/Radio-Calibration-BSL/tree/main/BSL) code involving the BSL method found including an example 'BSL.m'.

In [`/data_files`](https://github.com/HolgerBovbjerg/Radio-Calibration-BSL/tree/main/data_files) some simulated data files used in the scripts are found. These are mainly used in places were it would require extensive computation to generate the data.

In [`/misc`](https://github.com/HolgerBovbjerg/Radio-Calibration-BSL/tree/main/misc) a function for testing whether or not a propsed parameter is inside the prior range is found.

In [`/model_simulkation`](https://github.com/HolgerBovbjerg/Radio-Calibration-BSL/tree/main/model_simulation) a function for simulation of the Turin model is found including an example ´turin_simulation_theoretical_compare.m´. Both a gpu accelerated version and a non-gpu version is found.

In [`/summary_statistics`](https://github.com/HolgerBovbjerg/Radio-Calibration-BSL/tree/main/summary_statistics) the function for calculating summary statistics is found as well as an example of how the summary statistics change with the parameters.  

In [`/synthetic_likelihood`](https://github.com/HolgerBovbjerg/Radio-Calibration-BSL/tree/main/synthetic_likelihood) the function for calculating synethetic likelihood is found along with an example where the likelihood is calculated varying one of the parameters.  
