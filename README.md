[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5838392.svg)](https://doi.org/10.5281/zenodo.5838392)

# vahmani-etal_2021_AH

**Anthropogenic heating of the urban environment: an investigation of feedback dynamics between urban micro-climate and decomposed anthropogenic heating from buildings**

Pouya Vahmani<sup>1\*</sup>, Xuan Luo<sup>1</sup>, Tianzhen Hong<sup>1</sup>, Andrew Jones<sup>1</sup>

<sup>1 </sup> Lawrence Berkeley National Laboratory, Berkeley, CA 94720, USA

\* corresponding author:  pvahmani@lbl.gov

## Abstract
Cities consume 2/3 of global energy and consequently release a large amount of anthropogenic heat into urban environments, which are already vulnerable to extreme heat risk due to the compounding effects of urban heat island and the warming climate.  In this study, we use detailed process-based building energy modeling of over 1.1 million buildings in Los Angeles along with a high-resolution urban micro-climate modeling framework to assess the implications of anthropogenic heating for urban micro-climate dynamics and the feedback process between them. We uniquely distinguish between two major components of anthropogenic heating from HVAC system rejections and exhaust/relief air from buildings. We show that the less-studied anthropogenic heating from building exhaust, compared to that from HVAC systems, is mostly a nocturnal phenomenon with more significant implications for local air temperature due to the lower and more stable planetary boundary layer at night. We demonstrate that anthropogenic heating from HVAC rejection and building exhaust not only have reverse diurnal profiles, they also exhibit offsetting behaviors under increasing outdoor temperatures. Our results show that a detailed understanding of the composition of anthropogenic heating, specific to an urban environment, is required to predict its’ diurnal dynamics and its’ response to a warming climate.

## Journal reference
TBD

## Code reference
1. Vahmani, Pouya, Rastogi, Deeksha, & Thurber, Travis. (2021). IMMM-SFA/wrf_historical: v1.0.0 (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.5748047

## Data reference

### Input data
1. Rastogi, D., Vahmani, P., & Jones, A. (2021). ERA5-based 12-km WRF CONUS Historical Climate Data. [Data set]. Globus. https://app.globus.org/file-manager?origin_id=c296b088-b769-11eb-afd8-e1e7a67e00c1&origin_path=%2F.
2. Vahmani, P. et al. (2021). Input data for Vahmani et al 2021 AH. [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5834929.

### Output data
3. Vahmani, P. et al. (2021). Output data for Vahmani et al 2021 AH part 1 of 5. [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5838306.

4. Vahmani, P. et al. (2021). Output data for Vahmani et al 2021 AH part 2 of 5. [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5839178.

5. Vahmani, P. et al. (2021). Output data for Vahmani et al 2021 AH part 3 of 5. [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5847487.

6. Vahmani, P. et al. (2021). Output data for Vahmani et al 2021 AH part 4 of 5. [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5855041.

7. Vahmani, P. et al. (2021). Output data for Vahmani et al 2021 AH part 5 of 5. [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5855061.

## Contributing modeling software
| Model | Version | Repository Link | DOI |
|-------|---------|-----------------|-----|
| WRF | 4.0.1 | https://github.com/wrf-model/WRF/releases/tag/v4.0.1 | NA |
| WRF | 4.2.1 | https://github.com/wrf-model/WRF/releases/tag/v4.2.1 | NA |

## Reproduce my experiment
1. Follow the steps at IMMM-SFA/wrf_historical (https://doi.org/10.5281/zenodo.5748047) to produce historical data to be used as boundary conditions. The full data set is available at __[1]__.
2. Use the NDOWN package which is part of WRFv4.2.1 to downscale the results of step 1 over Los Angeles:
   * Follow the user guide for NDOWN: https://www2.mmm.ucar.edu/wrf/OnLineTutorial/CASES/NestRuns/ndown.php
   * Use the `namelist.input` and `namelist.wps` files provided at __[2]__.
   * The resulting subset of data is available at __[2]__.
3. Run the WU1 simulations:
   * Run WRF for the three nested domains over Los Angeles.
   * Use the `namelist.input`, all the `.TBL` files, the `myoutfields.txt` file, and the initial conditions `wrfbdy_d01`, `wrfinput_d01`, and `wrfinput_d02` provided at __[3]__ and __[4]__.
   * The results are available at __[3]__, __[4]__, __[5]__, __[6]__, & __[7]__.
4. Run the WU2 simulations:
   * Repeat the WU1 simulations, but with the anthropogenic heating parameters available in this repository here [here](3.AH_LA/) substituted into WRF's `URBPARM.TBL` parameter file.

## Reproduce my figures
Use the scripts from the folder [4.Figures/1.scripts](4.Figures/1.scripts) to reproduce the figures used in this publication:

| Script Name | Description | How to Run |
| --- | --- | --- |
| `M01_Variable_Extraction_for_hpc_V2.m` | Script to convert NetCDF files to MATLAB files | `matlab -nodisplay -nosplash -nodesktop -r "run('M01_Variable_Extraction_for_hpc_V2.m');exit;"` |
| `M01_Spatial_WRFout_Matlab_v5.m` | Script to create spatial maps in Figures 5 and S5 | `matlab -nodisplay -nosplash -nodesktop -r "run('M01_Spatial_WRFout_Matlab_v5.m');exit;"` |
| `M02_Temporal_WRFout_Matlab_v4.m` | Script to create timeseries plots in Figures 4, 6, and S3 | `matlab -nodisplay -nosplash -nodesktop -r "run('M02_Temporal_WRFout_Matlab_v4.m');exit;"` |
| `M03_Scatter_WRFout_Matlab_v4.m` | Script to create scatter plots in Figures 4, 7, 9, and S4 | `matlab -nodisplay -nosplash -nodesktop -r "run('M03_Scatter_WRFout_Matlab_v4.m');exit;"` |
