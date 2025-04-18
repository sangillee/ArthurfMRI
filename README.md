# ArthurfMRI

Released Apr 16th 2025

Set of functions that Arthur wrote to perform fMRI GLM analysis as fast as possible on MATLAB.
Adding each of the following four directories to MATLAB path should be enough for installation.
Documentations are a bit sparse at the moment, so I don't recommend using these yet...

glm: Contains functions for simulating BOLD responses from fsl-type regressors and running regression on multiple voxels simultaneously.

fmriprep: Contains function for running glm while interfacing with fmriprep preprocessed data. Essentially fmriprepRegression is the only user-facing function, and others are just helpers.

IO: contains two functions, Niftiopen and Niftisave, which builds upon MATLAB's niftiread and niftiwrite by performing numerical scaling when opening nifti files and automatically finding common nifti templates when writing. Particularly useful for directly saving vector-formatted brain data without explicitly loading mask images.

permutation: contains two functions for performing permutation tests.

searchlight: contains two functions for performing searchlight MVPA analysis.