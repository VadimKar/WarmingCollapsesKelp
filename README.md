# WarmingCollapsesKelp
Data and code to reproduce manuscript analyses

README FOR CODE & DATA BASE FOR KARATAYEV ET AL 
"Warming collapses giant kelp forests by depleting ecosystem resilience"



Code and Data provided here are solely for reproducing analyses.
They are not intended to serve as supplemental methods or new datasets.
Please do not re-use or re-distribute this data. Original can be found datasets in cited, openly accessible sources. 

Code contained herein involves several scripts which were run in sequence.
Each script generates data used in subsequent scripts and stores it in the "ProcessedData" directory.
For convenience, these processed data are already provided here, such that users may run Script 1 and then Script 5.

Script list:

Script 1: KelpAllee CoreFunctions.R
Enters core functions to use throughout analyses (necessary to run all other scripts)

Script 2: KelpAlleeBuildData0.3 LITE.R
Assembles diver and satellite data at LTM reefs

Script 3: GLM_KelpDynamics_0.3 LITE.R
Assembles region-wide remote sensing data and runs analyses in Figure 1

Script 4: KelpAlleeODEfit_0.5 ClimGen LITE.R
Builds future climate projections

Script 5: KelpAlleeODEfit_0.65 LITE.R
Runs model fitting, analysis, and validation (Figures 2-4)

Script 6: Kelp Allee ManyBfns Appendix.R
This is a standalone script to reproduce Methods




