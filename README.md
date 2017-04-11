# CVIVanalyzer
Extracts the electrical characteristics of a semiconductor device from a set of capacitance vs. voltage and current vs. voltage measurements

CV-IV Analyze instructions:

1.) Open "Configure.py" file with IDLE or Enthought.

2.) Enter data file locations and configuration parameters 
(Note: If Enthought is used, enter the full path for the 
locations of data files and the location where the plots to be saved)

3.) Run the program (If IDLE is being used, press f5 on Windows. On
Linux machines, open the terminal and cd to directory where the
configuration file is. Then, execute

> python nameofconfigurationfile.py


4.) Fitting ranges must be within the voltage range in the data.

5.) To fine-tune the fits, look at the pop-up plots and determine 
the fitting ranges by moving the cursor over the figure. The position 
of the cursor is shown at the lower right corner of the figure.

6.) The module "CVIVAnalyze.py" should be in one of the directories 
listed in "sys.path". To find these directories, do following
by opening a python shell:

> import sys
> sys.path

This will list the directories in the python system path.

7.) Keep configuration files for each measurement by giving unique 
file names. After finishing the current analysis, keep the configuration 
file and make another file to do the next analysis. This will save a lot 
of time in case the analysis is needed to be revisited.
