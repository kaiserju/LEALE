This repository hosts LEALE, a tool designed to estimate the salt wedge intrusion length in estuarine and riverine systems.

To use this tool, the first step to be done is run ETICO tool (https://github.com/fabioviola-cmcc/EstuarIO_thalweg.git) to get the thalweg path.
Within this file, you should make sure you have model output (from where you want to get the salinity information) in netCDF format.

Inputs from the user are:
- The complete filename from ETICO (pathnameETICOfile);
- Pathname where is the model output (pathnamein);
- Pathname where the result from LEALE will be saved (pathnameout)

To run the script, the following packages are needed: os, netCDF4, pandas, pyresample, geopy, and datetime.
