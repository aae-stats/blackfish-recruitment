# VEFMAP Stage 7: responses of blackfish recruitment to flow

## Description

This repository contains code to analyse associations between recruitment of 
blackfish (*Gadopsis marmoratus* and *Gadopsis bispinosus*) and a set of flow
and water temperature variables. Recruitment was defined as survival to 
young-of-year and 1+ age classes, with age estimated from length data using
length thresholds (< 80 mm for young-of-year, 80-160 mm for 1+).


## Usage

The entire analysis workflow is contained in `main.R`, supported by a series
of helper functions in the `R` directory. Outputs are saved to an `outputs` 
directory, and all data is loaded from a database and cached in the `data`
directory. Data and outputs are not included here; outputs are provided
in the final VEFMAP report (Amtstaetter et al. 2024) and data are available
on request.


## Additional information

Queries can be directed to jian.yen [at] deeca.vic.gov.au.


## References 

Amtstaetter, F., Tonkin, Z., Yen, J., Lieschke, J., Sharley, J., Raymond, S., Koster, W., Harris, A., Crowther, D. and Ryan, S. (2024). Flow regime and water temperature influence recruitment success of Blackfish in south-eastern Australia. Arthur Rylah Institute for Environmental Research Unpublished Client Report. Department of Energy, Environment and Climate Action, Heidelberg, Victoria. 
