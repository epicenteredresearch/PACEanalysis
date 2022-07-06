# PACEanalysis

## UPDATES IN VERSION 0.1.8
- Updates to ensure compatible with updates to dependencies 
- Updating code for sesame detection p-values due to changes in sesame package
- Adding option to run site-specific analysis in parallel to speed up computation time, see dataAnalysis function for more details
- Updating analysis of the association between cell composition estimates and model variables

## UPDATES IN VERSION 0.1.7
- Now creates two "Table 1"s - one for all samples, and one complete case version
- Adding cell composition to Table 1 (if Omega is provided)
- Updated the naming of the Table 1 file so that multiple versions can be open at the same time

## UPDATES IN VERSION 0.1.6
- Offers two options to correct for the detection p-values; the default is the SeSAMe approach (same as before)
- Offers investigators the option to choose their detection p-value cut-off
- Offers the option to remove probes with less than a certain number of beads (recommended minimum number is 3)
- Offers the option to remove probes with zero intensity values
- Offers the option to choose whether to trim or winsorize the outliers. If winsorizing, the investigator can choose the percent to winsorize on each end. If trimming, removes extreme outliers based on gaps that must be at least 3*IQR; the cutoff for the number of outliers in a group is either a maximum of 5 or 0.0025 of the total number of samples (whichever is larger). 
- Corrected the issue with -Inf and Inf estimates of contamination
- Made sure that the package is compatible with the recent R update (v4.1.0)
- Masking poor probe intensities has been moved into a new function so that investigators can change their mind without rerunning the full preprocessing function
- Dealing with outliers has been moved into a new function so that investigators can change their mind without rerunning the full preprocessing function
