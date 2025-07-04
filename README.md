This is a code repository for the manuscript “Modeling age patterns of under-5 mortality: a comprehensive model life table approach for low- and middle-income countries”, which is currently under review.

The MATLAB file Model_B.m produces the figures and tables included in the manuscript, using one input data file: BoOTStrAp_DHS.mat. This data file is not part of this repository but can be generated after running the Stata do-file _DHS_processing.do and the MATLAB file Model_B_BoOTStrApPInGDHS.m in that order (both files within this repository). 

These routines require FLat data files: XXIRYYFL.dta (individual recode for all women in reproductive ages) and XXBRYYFL.dta (birth recode consisting of full birth histories) by the DHS Program; where XX indicates the country prefix (as described in _codesDHS.dta), and YY is the corresponding survey suffix of all Demographic and Health Surveys available at the time of writing this manuscript (as described in lISt.csv). These FLat files are available—upon request—at: https://dhsprogram.com

The do-file and m-files run automatically from top to bottom, but the user may need to adjust the file paths for reading the data and saving the outputs. The m-files may require some nested functions (also within this repository) to do some part of the analysis. Additional nested m-functions to produce tables and for under-5 mortality modelling are available at: https://github.com/Romero-Prieto/RaMMPS_U5M (also within this GitHub account).   
# Model_B
