# README

This repository contains the code and data for the research paper by Lin et.al.

## Directory Structure and File Descriptions

- **codes/**  
  Contains the Random forest models used in this study for predicting the parameter N.
  Also contains the code for implementing baseflow separation using the SMM method.

- **Mei et.al data/**  
  - **SO_FSpl-SMM-EC/**  
    Contains the runoff data for the basins used in this study.  
  - **ResTab.mat**  
    Contains the optimal parameter values (N) for the SMM for each basin.

- **RF2.pkl, RF3.pkl, RF4.pkl**  
  Random forest models used in this study for predicting the parameter N.

- **run_Nprediction.py**  
  Python script for predicting the parameter N using the provided random forest models.

- **stationInfo.xlsx**  
  Input file required by `run_Nprediction.py`. It contains the attribute features for each basin.

- **predicted_N.xlsx**  
  Output file generated after running `run_Nprediction.py`, which contains the predicted parameter N values for each basin.

- **run_SMM.m**  
  MATLAB script for performing baseflow separation using the SMM method. This script requires the `predicted_N.xlsx` file along with the runoff data.

- **BF_result/**  
  Directory that stores the baseflow separation results produced after running `run_SMM.m`.

## Usage

1. **Parameter Prediction:**  
   Run `run_Nprediction.py` using the `stationInfo.xlsx` file to predict the optimal parameter N values. The results will be saved in `predicted_N.xlsx`.

2. **Baseflow Separation:**  
   Execute `run_SMM.m` in MATLAB, ensuring that `predicted_N.xlsx` and the required runoff data from `Mei et.al data/SO_FSpl-SMM-EC/` are available. The output baseflow separation results will be stored in the `BF_result/` directory.

