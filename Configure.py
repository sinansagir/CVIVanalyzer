#!/usr/bin/python

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""     CONFIGURATION PARAMETERS    """"""""""""""""""""""""
""""""""""""""""""""""""          VERSION : 3.17         """"""""""""""""""""""""
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

CVFileLocation = "/Users/sinansagir/Research/Silicon_Research/CVIVAnalysis/CVIVAnalyze3.17/Data/MCZ200Y_06_DiodeL_5_7_5_2013_12_54_22_ PM.txt" %locals()
IVFileLocation = "/Users/sinansagir/Research/Silicon_Research/CVIVAnalysis/CVIVAnalyze3.17/Data/MCZ200Y_06_DiodeL_5_7_5_2013_12_45_09_ PM.txt" %locals()
PlotSavingLocation = "/Users/sinansagir/Research/Silicon_Research/CVIVAnalysis/CVIVAnalyze3.17/Plot/" %locals()

d = 200.0 * 1e-6 # thickness of sensor in SI measured before irradiation (if non-irradiated, it will be extracted from measurement via end capacitance)
A = 25.0 * 1e-6 # active area of sensor in SI
Irradiation = 0 # 1 for irradiated samples, 0 for non-irradiated samples
SensorType = 1 # 1 for P-type, 0 for N-type
ScaleTemp = -20.0 # Current scaling temperature
FlatFittingMethod1 = 1 # 0 for y=ax+b and 1 for y=b (forced to be horizontal line) ## for log-log scaled plot
FlatFittingMethod2 = 1 # 0 for y=ax+b and 1 for y=b (forced to be horizontal line) ## for 1/C^2 plot
CapManualf1 = -12.85E-12 # F (>0 if manual capacitance is wanted for flat fit and <0 if manual capacitance is not wanted)
CapManualf1Err = 10E-15 # F
CapManualf2 = -12.77E-12 # F (>0 if manual capacitance is wanted for flat fit and <0 if manual capacitance is not wanted)
CapManualf2Err = 10E-15 # F

slopef1min = 25
slopef1max = 40
flatf1min = 200
flatf1max = 400

slopef2min = 25
slopef2max = 45
flatf2min = 200
flatf2max = 400

"""If different ranges are wanted for log-log and 1/C^2 fittings, use following
parameters in addition to the parameters above (in this case, above parameters
set the ranges for log-log scaled fittings). If the same ranges are wanted for
both log-log and 1/C^2 fittings, then set slopef1min_1C < 0 and use the above
parameters to set the fitting ranges!"""

slopef1min_1C = -170
slopef1max_1C = 220
flatf1min_1C = 300
flatf1max_1C = 600

slopef2min_1C = 150
slopef2max_1C = 200
flatf2min_1C = 300
flatf2max_1C = 600

Vmin = 2 
Vmax = 1000

CapMode = 0 # 0 for parallel mode, 1 for serial mode
CapOffset = 0 # Capacitance offset if any
VdepCoeff = 1.2 # Vdep coefficient for extraction of leakage current and capacitance above full-depletion (for HPK, it is 1.2)
PlotDopingProfile = 0 # 0 if plot not wanted, 1 if plot wanted
PlotDoubleDerivative = 0 # 0 if plot not wanted, 1 if plot wanted
SavePopUpPlots = 0 # 0 if pop-up plots not wanted to be saved, 1 if wanted
LCRAccuracy = 0.05 # %
OpenMeasError = 3e-14 # F
CurrentAccuracy = 0.1 # %

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Sensor          Thickness[um]  Active Area[mm2]  Capacitance[pF]
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
CA0126              284.7           25.0              9.252 
CA0127              285.8           25.0              9.216
CA0133              255.4           25.0              10.315

CD1731              282.4           25.0              9.328
CD1737              280.2           25.0              9.400
CD1738              289.8           25.0              9.091

8556-02-37-3        300.0           6.25              2.880
8556-02-37-4        300.0           6.25              2.877
8556-02-39-1        300.0           6.25              2.876
8556-02-39-2        300.0           6.25              2.877

261636-10-23-1      155.1           6.25              4.246
261636-10-23-3      154.2           6.25              4.270
261636-10-23-4      155.4           6.25              4.238
261636-10-29-1      153.6           6.25              4.287

6336-03-41          49.11           25.0              53.642
6336-03-42          49.51           25.0              53.209
6336-03-43          49.85           25.0              52.836
6336-03-44          51.03           25.0              51.616

MCZ200Y_06_DiodeL_5 200.0           25.0
FTH200Y_04_DiodeL_9 200.0           25.0
MCZ200N_06_DiodeL_9 200.0           25.0
FTH200N_25_DiodeL_8 200.0           25.0
MCZ200Y_04_DiodeL_3 200.0           25.0
FTH200Y_03_DiodeL_9 200.0           25.0
MCZ200N_11_DiodeL_9 200.0           25.0
FTH200N_25_DiodeL_3 200.0           25.0
MCZ200Y_07_DiodeL_9 200.0           25.0
FTH200Y_04_DiodeL_8 200.0           25.0
MCZ200N_09_DiodeL_8 200.0           25.0
FTH200N_25_DiodeL_2 200.0           25.0
MCZ200Y_05_DiodeL_2 200.0           25.0
FTH200Y_04_DiodeL_3 200.0           25.0
MCZ200N_11_DiodeL_8 200.0           25.0
FTH200N_24_DiodeL_8 200.0           25.0
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

from CVIVAnalyze import *

ConfigurationParameters(CVFileLocation, IVFileLocation, PlotSavingLocation, d, A, Irradiation, SensorType, ScaleTemp, FlatFittingMethod1, FlatFittingMethod2, CapManualf1, CapManualf1Err, CapManualf2, CapManualf2Err, slopef1min, slopef1max, flatf1min, flatf1max, slopef2min, slopef2max, flatf2min, flatf2max, slopef1min_1C, slopef1max_1C, flatf1min_1C, flatf1max_1C, slopef2min_1C, slopef2max_1C, flatf2min_1C, flatf2max_1C, Vmin, Vmax, CapMode, CapOffset, VdepCoeff, PlotDopingProfile, PlotDoubleDerivative, SavePopUpPlots, LCRAccuracy,OpenMeasError, CurrentAccuracy)

