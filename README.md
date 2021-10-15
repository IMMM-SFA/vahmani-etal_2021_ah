_your zenodo badge here_

# vahmani-etal_2021_tbd

**Anthropogenic heating of urban environment: an investigation of feedback dynamics between urban micro-climate and components of anthropogenic heating from buildings**

Pouya Vahmani<sup>1\*</sup>, Xuan Luo<sup>1</sup>, Tianzhen Hong<sup>1</sup>, Andrew Jones<sup>1</sup>

<sup>1 </sup> Lawrence Berkeley National Laboratory, Berkeley, CA 94720, USA

\* corresponding author:  pvahmani@lbl.gov

## Abstract
TBD

## Journal reference
TBD

## Code reference
Currently here for workflow code but needs to be confirmed:  https://github.com/IMMM-SFA/phase_2_climate

## Data reference

### Input data
WRF CONUS sims on globus, boundary conditions come from TGW-WRF
Just historic data.

### Output data
WRF outputs for ~400m resolution.  Currently on Cori.

## Contributing modeling software
| Model | Version | Repository Link | DOI |
|-------|---------|-----------------|-----|
| WRF | 4.2.1 | https://github.com/wrf-model/WRF/releases/tag/v4.2.1 | NA |

## Reproduce my experiment
Currently here for workflow code but needs to be confirmed:  https://github.com/IMMM-SFA/phase_2_climate

## Reproduce my figures
Use the scripts found in the `figures` directory to reproduce the figures used in this publication.

| Script Name | Description | How to Run |
| --- | --- | --- |
| `generate_figures.py` | Script to generate my figures | `python3 generate_figures.py -i /path/to/inputs -o /path/to/outuptdir` |
